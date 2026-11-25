# Individuelle referencer og pooled BAM

Disse selvstændige scripts er **kandidater til forskning og validering**, ikke
klinisk validerede. De ændrer ikke eksisterende RefSample eller v6.4.
Kontrollér prøvernes identitet, behandling, assay, reference, QC,
slægtskab og genotype før brug. `qc=pass` og brugerens valg er ikke
uafhængig kvalitetskontrol.

Krav: `Rscript` med base-R, `samtools` og `regtools` på PATH. Der bruges
indekserede coordinate-sorterede BAM'er med samme sekvensdictionary i
samme rækkefølge. Ens contignavne og -længder beviser ikke, at alignment,
genomsekvens eller øvrige protokoller er identiske. Brug kun individuelt
godkendte, sammenlignelige BAM'er fra samme væv og protokol.

## 1. Lav en pool fra en BAM-liste

Opret to tekstfiler, **én for hver betingelse**, med én absolut BAM-sti pr.
linje. Her er stierne kun eksempler og skal erstattes med faktiske filer:

```text
/drift/RNA/reference_inputs/patient1_nt.bam
/drift/RNA/reference_inputs/patient2_nt.bam
```

Tomme linjer og linjer med `#` først ignoreres. Samme BAM-sti må ikke
optræde to gange. Tjek manuelt at hver liste kun indeholder den angivne
betingelse og ikke patienten, der skal sammenlignes; scriptet kan **ikke**
udlede behandling eller person fra BAM-navnet. Inkludér kun én relevant
prøve pr. person. De to filer må ikke være en samlet liste af nt og puro.

Eksempel med separate outputmapper, som ikke må eksistere i forvejen:

```bash
Rscript /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/scripts/reference_pool.R \
  --list /drift/RNA/reference_inputs/nt_bams.txt \
  --condition nt \
  --output-dir /drift/RNA/TotalRNA_Splice/RefSample/ClinV1_nt_v1 \
  --target-reads 20000000 --seed 100 --threads 4

Rscript /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/scripts/reference_pool.R \
  --list /drift/RNA/reference_inputs/puro_bams.txt \
  --condition puro \
  --output-dir /drift/RNA/TotalRNA_Splice/RefSample/ClinV1_puro_v1 \
  --target-reads 20000000 --seed 100 --threads 4
```

`--target-reads` er et **omtrentligt mål for mappede alignments pr.
input-BAM**, ikke et garanteret antal reads i output. Faktoren beregnes ud
fra `samtools idxstats`; sekundære/supplementære alignments kan påvirke
tallet, og subsampling er stokastisk og udføres på templates. BAM'er med
færre mappede alignments end målet bruges uden subsampling; unmapped
alignments er ikke særskilt filtreret. Dette er ikke normalisering af
genekspression. Standard er 20 millioner, seed 100 og 2 tråde.
Ønskes alle reads anvendt, vælg et mål større end antallet af mappede
alignments i hver BAM. Der køres én fil ad gangen; større `--threads` er
primært til samtools-komprimering.

Resultater: `reference_nt.bam` eller `reference_puro.bam` og BAI,
`inputs.tsv` med input, faktor og filmetadata, sekvensdictionary,
samtools-version, logs og `STATUS.txt`. Midlertidige subsamples bliver
liggende i `subsamples/`; planlæg diskplads. Ved fejl står status som
`INCOMPLETE`, og mappen overskrives ikke. Undersøg logs og brug en ny mappe.
`quickcheck` læser ikke hele filen og garanterer ikke mod intern korruption.
En pool er ét visuelt spor, ikke uafhængige kontroller; fjern/udelad den
aktuelle person ved opbygning af den pool, der bruges til vedkommende.

ClinV1's `-refs` anvender **alle** BAM-filer direkte i den angivne mappe
til **begge** patientbetingelser. Det kan derfor ikke automatisk matche
disse separate pools korrekt. Åbn dem separat i IGV eller kør uden `-refs`;
angiv ikke én pool som betingelsesmatchet kontrol for begge prøver.

## 2. Individuel junction-referencebank: én beregning pr. prøve

Dette er en selvstændig referencebank, **ikke** en automatisk del af
ClinV1-kørslen. Den kører `regtools junctions extract -s XS -a 8 -m 50
-M 500000` én gang på hele hver input-BAM og gemmer individuelle tællinger
opdelt på kromosom. Det er en tung første beregning, ikke blot et
genudsnit. Ved nye prøver beregnes kun de nye BAM'er. Sammenligning
læser kun det relevante kromosom fra de valgte prøver.

Den måler **kun junctions**: ingen vurdering af intronretention,
lokal dækning, allelspecifik ekspression, normaliseret splejsningsandel,
p-værdi eller funktionel effekt. Den er heller ikke en implementering af
FRASER/OUTRIDER og er ikke nødvendigvis direkte kompatibel med deres input.
XS-tags og regtools-filtre har samme begrænsninger som i ClinV1. Værdierne
`.`, `*` og `?` bevares hver for sig som ukendt/uafklaret strand og matches
ikke med `+` eller `-`; kontrollér strand og bibliotek før fortolkning.

### Opret banken

`--reference-id` er en brugerangivet identifikator for den fastlagte
alignment/referenceprotokol, ikke et automatisk bevis for referenceversion:

```bash
Rscript /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/scripts/reference_bank.R \
  init --bank /drift/RNA/TotalRNA_Splice/RefSample/ClinV1_bank_v1 \
  --reference-id GRCh38_totalRNA_v1
```

Banken fastlåser regtools-binærens MD5 og parametrene. Hvis værktøjet
skiftes, må prøver ikke tilføjes til samme bank: opret en ny version og
genberegn/validér. BAM-contigernes navne, længder og rækkefølge tjekkes
mod første prøve; andre forhold kræver manuel dokumentation. Gem bank,
manifest og BAM-stier under adgangskontrol.

### Tilføj prøver

Opret en **tabulatorsepareret** manifestfil med præcis disse kolonner:

```text
sample_id	person_id	condition	bam	qc
P01_nt	P01	nt	/drift/RNA/reference_inputs/P01_nt.bam	pass
P01_puro	P01	puro	/drift/RNA/reference_inputs/P01_puro.bam	pass
P02_nt	P02	nt	/drift/RNA/reference_inputs/P02_nt.bam	pass
```

Erstat eksemplerne; `person_id` skal være den samme for begge konditioner
fra samme person. `sample_id` og `person_id` skal være pseudonyme,
entydige simple ID'er med bogstaver, tal, `_`, `-` eller `.`. `condition`
skal være `nt` eller `puro`; `qc` skal være `pass`. BAM-stier skal være
absolutte og have BAI/CSI. Én person kan højst have én prøve per
betingelse, og samme BAM må ikke genbruges som en anden prøve.

```bash
Rscript /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/scripts/reference_bank.R \
  add --bank /drift/RNA/TotalRNA_Splice/RefSample/ClinV1_bank_v1 \
  --manifest /drift/RNA/reference_inputs/manifest.tsv
```

Et nyt manifest kan indeholde kun de nye prøver eller alle hidtil
registrerede prøver. Uændrede, tidligere tilføjede prøver springes over;
ændret ID, betingelse, BAM-sti, størrelse eller ændringstid afvises.
Filstørrelse og tidspunkt er ikke en kryptografisk integritetskontrol.
Prøven publiceres først efter vellykket tælling. Hvis flere prøver tilføjes
i samme kørsel, kan tidligere prøver være publiceret, selvom en senere
fejler; ret fejlen og genkør manifestet. Fejlede staging-filer bevares til
gennemgang. Kun én `add`-proces må skrive ad gangen; ved et brat afbrudt
job kan `.add.lock` blive liggende. Fjern kun låsen manuelt, når det er
bekræftet, at ingen skriver til banken.

### Sammenlign en patient mod udvalgte personer

Lav en tekstfil med præcis de **sample_id'er**, der efter manuel vurdering
er egnede referencer for det aktuelle gen og den aktuelle betingelse,
ét ID pr. linje. Brug ikke den aktuelle patient, en anden prøve fra samme
person eller uegnede/relaterede personer. Genotype- og fænotypefiltrering
sker **ikke automatisk**. Eksempel:

```text
P02_nt
P03_nt
```

Brug patientens rå junction-BED fra en fuldført ClinV1-kørsel og angiv
den undersøgte region som **1-baseret inklusivt** `chr:start-end`.
Fil- og regionsnavne nedenfor er pladsholdere:

```bash
Rscript /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/scripts/reference_bank.R \
  compare \
  --bank /drift/RNA/TotalRNA_Splice/RefSample/ClinV1_bank_v1 \
  --patient P01 --condition nt \
  --junctions /drift/RNA/path/to/P01_nt.bed \
  --region chr22:START-END \
  --controls /drift/RNA/reference_inputs/CHEK2_nt_controls.txt \
  --output-dir /drift/RNA/TotalRNA_Splice/TestRuns/P01_CHEK2_nt_refs
```

Gentag med `puro`, puro-BED og en særskilt puro-kontrolliste.

Gamle `minus`/`plus`-værdier accepteres fortsat som alias for `nt`/`puro`.
Eksisterende bankmetadata læses uden migrering; nye metadata bruger `nt`/`puro`.
Ved pool-kommando med et gammelt alias bevares det gamle outputfilnavn,
fx `reference_minus.bam`; med `nt` hedder det `reference_nt.bam`.
Sørg for, at patientens BED er lavet med samme regtools-version, parametre,
reference og regionsafgrænsning; dette kan scriptet ikke verificere ud
fra BED-filen. Patient-BED fra et afgrænset BAM-udsnit kan mangle reads
ved regionens grænser, som findes i referencebankens fulde BAM.

`comparison.tsv` har én række for hver junction, der overlapper regionen
og er observeret i patienten **eller** en valgt kontrol, ganget med valgte
kontroller. Junctioner matches med eksakte 0-baserede, halvåbne
introngrænser og strand. `not_observed_coverage_unassessed` og `NA` betyder
ikke, at hændelsen er biologisk fraværende. Rå counts er ikke direkte
normaliseret mellem prøver. Patientens BED kan være tom, og et tomt
resultat betyder heller ikke, at området kan vurderes.

`controls_used.tsv`, `patient_junctions.tsv`, `control_inventory.tsv` og
`STATUS.txt` dokumenterer sammenligningen, herunder filchecksums og
reference-ID. Bevar resultater og bankversion uændret; en senere tilføjelse
af prøver giver ikke automatisk samme kontroludvalg. Brug en ny outputmappe
ved genkørsel. Kun registrerede kromosomfiler med korrekt checksum læses;
fejlede filer må ikke fortolkes som nul.

## Tests

```bash
Rscript /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/tests/test_reference_pool.R
Rscript /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/tests/test_reference_bank.R
bash /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/tests/test_reference_workflows.sh
```

De to første bruger simulerede eksterne værktøjer. Den sidste opretter
syntetiske BAM-filer og bruger de rigtige lokale samtools og regtools.
Ingen af testene er biologisk eller klinisk validering.
