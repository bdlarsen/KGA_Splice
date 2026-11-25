# ClinV1: kandidat til målrettet RNA-spliceanalyse

ClinV1 er en separat, reduceret version. **Den er ikke klinisk valideret eller
godkendt og må ikke bruges alene til patientbeslutninger.** Kør ikke ClinV1
samtidig på samme outputmappe som en anden kørsel. v6.4 ændres ikke.

## Kom godt i gang

På Docker-værten kan den faste ClinV1-opsætning startes med:

```bash
bash /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/run_clinv1.sh
```

Scriptet spørger om variant, gen, transkript, prøve-ID'er, genotype og
BAM-mappe. `--help` viser alternativet med parametre og `--yes`. De faste
referencer, SAI, image-ID og outputplacering er defineret i scriptet; et
ændret image-ID skal kontrolleres og valideres før scriptet opdateres.
Resultater skrives som standard direkte i `/drift/RNA/drift/ID_GEN_unik/`,
herunder `tmp`. En anden eksisterende, skrivbar outputmappe kan vælges
interaktivt eller med `--output /absolut/sti`. Tryk Enter for standardmappen.
Værtens to opdaterede R-scripts mountes read-only oven på image-kopierne;
deres SHA-256 og image-ID gemmes i `tmp/launcher_provenance.txt`.
Tidligere resultater overskrives ikke. Containeren kører uden netværk,
input mountes read-only og stierne er de samme på vært og i IGV-session.
Kontrollér plads på `/` og `/drift`: ikke alle containerfiler ligger i `tmp`.

### IGV på værten og i VSwarehouse

Hver ny kørsel opretter to sessioner i projektets `IGV`-mappe:

* `paired_IGV.xml` til IGV, der ser filer under værtens `/drift`.
* `paired_IGV_vswarehouse.xml` til IGV i VSwarehouse, der ser samme drev
  under `/home/ghuser/Workspace/drift`.

Den anden session oversætter ressource- og trackstier fra `/drift/...` til
`/home/ghuser/Workspace/drift/...`. Den ændrer ikke selve BAM-filerne.
Kontrollér i VSwarehouse, at både BAM, BAM-indeks og projektets output kan
læses på de oversatte stier. Et eksisterende indeks på analyseværten er ikke
tilstrækkeligt, hvis VSwarehouse ikke kan se det. Stier uden for `/drift`
oversættes ikke: pipelinen giver en advarsel, og disse tracks kan kræve en
anden mount eller manuel tilpasning. IGV skal også have adgang til hg38.

En tidligere kørsel kan få en VSwarehouse-session uden at genkøre analysen.
På Docker-værten angives først den eksisterende projektmappe; kommandoen
overskriver ikke en eksisterende VSwarehouse-fil:

```bash
PROJECT=/drift/RNA/drift/UDSKIFT_MED_PROJEKTMAPPE
SCRIPTS=/drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/scripts

docker run --rm --network none \
  --mount "type=bind,src=$PROJECT,dst=/project" \
  --mount "type=bind,src=$SCRIPTS,dst=/clinv1-scripts,readonly" \
  kga-splice:clinv1-minimal-candidate \
  Rscript /clinv1-scripts/make_vswarehouse_session.R \
  /project/IGV/paired_IGV.xml /project/IGV/paired_IGV_vswarehouse.xml
```

Åbn derefter `paired_IGV_vswarehouse.xml` i VSwarehouse. Verificér at alle
tracks indlæses; oprettelsen af XML-filen kan ikke kontrollere adgang i det
andet miljø. Til eventuelt IGV-snapshot bruger analysen fortsat den normale
session.

For nye kørsler sættes samme IGV *Junction Coverage Min* på alle patient-
og kontrolspor, både regions- og original-BAM. Tærsklen er 1 % af medianen
af de positive regtools-junctiontællinger i den udtrukne region (genet
+/- 10 kb) for første patientprøve (-puro), afrundet op og begrænset til
3–20 reads. Ved enkeltprøver med ukendt
behandling bruges denne prøve. Hvis den ikke har positive junctions,
bruges IGV's egen standard, og pipelinen giver en advarsel. Den valgte
tærskel eller manglende beregningsgrundlag dokumenteres i
`IGV/junction_threshold.txt`.

Dette er et eksplorativt visningsfilter, ikke en valideret analysegrænse.
Det ændrer ikke junctiontabellen; svagt understøttede junctions kan skjules.
Sænk grænsen manuelt i IGV, f.eks. til 1, når relevante svage junctions
undersøges. Regtools- og IGV-tællinger kan afvige. Konvertering af en
eksisterende session til VSwarehouse ændrer ikke dens tærskel.

Kør `Rscript scripts/ClinV1.R --help` (eller uden argumenter) for en kort
oversigt. `-variant` kan være den eksisterende `CHR-POS-REF-ALT`-værdi
eller stien til en `.vcf`-/`.vcf.gz`-fil med præcis én biallelisk variant.
Flere records, multialleliske og symbolske alleler afvises; lav ét separat
run pr. variant. VCF-samplekolonner og genotype anvendes ikke; `-gt` skal
angives separat, hvis patientens genotype er kendt. FASTA kontrollerer REF.
VCF'ens assembly verificeres ikke automatisk; input skal være GRCh38.
`variant_input.txt` og `STATUS.txt` angiver input og genomisk variant.
HGVS-input er endnu ikke understøttet.

Normalt opretter `-p` en ny projektmappe under den angivne mappe. Alternativt
kan `-project-dir` anvise en allerede reserveret mappe, som kun må indeholde
en valgfri `tmp`-undermappe; output skrives direkte i projektmappen. Angiv
præcis én af `-p` og `-project-dir`, og brug altid en ny mappe pr. kørsel.

`-s ID` kører én patientprøve med **ukendt behandling**; rapporten viser
ingen nt/puro-sammenligning. `-s ID_NT,ID_PURO` bevarer den eksisterende
parkørsel: første prøve uden og anden med puromycin. Behandling og identitet
skal verificeres uafhængigt. BAM-filer søges direkte i `-bam`-mappen med
ID som filnavnsstart, eventuelt med præfikset `G150-`; ID skal efterfølges
af `.bam`, `_` eller `-`. Flere matches afvises, og undermapper søges ikke.
Ved batchkørsel kan `scripts/batch_preflight.R` kontrollere outputmappe,
SAI-modul og patient-/reference-BAM-indekser før analysen. Scriptet læser
linjer `ID|variant|gen|transkript|BAM-mappe|single` eller `...|paired`
fra standard input. Det verificerer ikke prøveidentitet eller BAM-indhold.

### Valgfrit IGV-snapshot

Tilføj `-igv-snapshot /absolut/sti/til/igv.sh` for at forsøge at generere
`Report/igv_snapshot.svg` i et vindue på +/- 1000 bp omkring varianten.
Kræver IGV Desktop med `igv.sh`, `xvfb-run`, Java og en lokal indekseret FASTA.
Brug en netværksisoleret container (`--network none`), så hverken patientdata
eller varianter kan sendes over netværket. ClinV1 håndhæver ikke selv denne
isolering. IGV kan have yderligere lokale referencekrav; installations- og
runtimekompatibilitet skal testes med syntetiske data først.

Funktionen er deaktiveret som standard. Fejl stopper ikke analysen, men
beskrives i konsollen, `STATUS.txt` og rapporten; se også
`Logs/igv_snapshot.out` og `Logs/igv_snapshot.err`. SVG-kontrol bekræfter
kun filformat, ikke at alle tracks er gengivet. Kontrollér billedet visuelt
mod IGV-sessionen; tracks kan mangle eller klippes ved mange prøver. Sessionens
alignment-tracks er skjult som standard; billedet er en statisk visning,
ikke en klinisk konklusion. Stier til IGV, FASTA og projekt må ikke indeholde
whitespace ved automatisk eksport. Eksisterende rapporter ændres ikke.

* [Enkel guide: start container og kør ClinV1](../DOCKER_USAGE.md#enkel-guide-clinv1-med-docker)
* [Docker: build, mounts, tests og kørsel](docs/DOCKER.md)
* [Separat minimalt ClinV1-image og testplan](docs/MINIMAL_DOCKER.md)
* [Git og GitHub: trin-for-trin](docs/GITHUB.md)
* [Individuel referencebank og BAM-pools](docs/REFERENCEBANK.md)
* [Kohorteplan og validering](docs/KOHORTEPLAN.md)

## Hvad køres?

1. Entydig identifikation af én patient-BAM med ukendt behandling eller to
   forskellige BAM-filer: nt uden puromycin, puro med puromycin. BAM-indeks,
   FASTA, GTF, transkript og REF-allel tjekkes.
2. Én fælles lokal VEP- og SpliceAI-annotering. Rå VCF bevares; manglende
   annotation er **ikke** en negativ prædiktion. VEP-output kan indeholde
   flere transkripter og skal fortolkes for det valgte transkript.
3. BAM-udsnit for genet +/- 10 kb via `samtools view` med BAM-indeks.
   `samtools quickcheck` kontrollerer inputfilens header og EOF uden en
   fuld scanning; den kan ikke opdage alle interne fejl. Regionsudtrækket
   fejler ved problemer med læsning af den anvendte region eller indeks.
   Øvrige dele af BAM-filen kontrolleres ikke for korruption.
   Kromosomnavn og -længde i BAM sammenholdes med FASTA; dette alene beviser
   ikke, at hele BAM'en og FASTA er fra samme referenceversion. GTF-filen
   gennemlæses for hvert run, hvilket kan tage tid på et netværksdrev.
4. Junctions fra hver prøve og eventuelle reference-BAM'er optælles med
   regtools. Eksakte introngrænser og strand sammenlignes med GTF-exons i
   det valgte transkript og øvrige annoterede transkripter i genet.
5. Variantpositionens dækning, strandfiltrerede bedgraphs og en fælles
   IGV-session med nt før puro og de originale BAM'er nederst.
6. Deskriptiv alleltælling ved en angivet heterozygot SNP og en samlet
   HTML-rapport med VEP, SpliceAI, dæknings-/junctionstatistik og figurer.

Det produceres **ikke**: automatisk vurdering af observeret pseudoexon,
intronretention, frameshift eller NMD-effekt, automatisk MAE-klassifikation,
signifikans, klinisk klassifikation eller valideret klinisk konklusion. Brug
IGV-læsninger og lokale protokoller til manuel vurdering; rå junctions
er ikke i sig selv en komplet analyse af disse effekter.

## Krav

Linux/Docker-miljø med `Rscript` og R-pakken `xml2`, `samtools`, `sambamba`
`bedtools`, `bash` og `regtools` på PATH samt Conda-miljøerne `vep_env` og `spliceai_env`.
Der kræves desuden `sashimipy` fra pakken `sashimi.py` i Conda-miljøet
`KGA_Splice`, som i v6.4. Tilgængelighed kontrolleres før annotering.
For `-gt het` på en SNP kræves også `gatk` og en FASTA-sekvensdictionary
som `genome.dict` ved siden af `genome.fa`.
Lokalt VEP RefSeq-cache, GRCh38 FASTA med `.fai`, GTF og coordinate-sorterede
BAM-filer med `.bai` eller `.csi` er nødvendige. De eksisterende standardstier
er de samme som i v6.4; de kan erstattes med `-g`, `-gtf` og `-vep`.
Kontrollér at regtools-versionen understøtter `-s XS`. SpliceAI's indbyggede
`grch38`-annotation og FASTA/BAM/GTF skal være indbyrdes kompatible.
ClinV1 kontrollerer ved start, at kommandoerne findes på PATH. Hvis
`regtools` mangler, skal det installeres eller gøres tilgængeligt i
containerens aktive miljø; se fejlmeddelelsen og kontrollér med
`command -v regtools`. Denne kontrol beviser ikke, at værktøjet kan køre
eller har en kompatibel version.

## Testkommando i containeren

Kør først den syntetiske test uden patientdata:

```bash
Rscript /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/tests/test_core.R
Rscript /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/tests/test_report.R
Rscript /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/tests/test_sai.R
Rscript /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/tests/test_report_pipeline.R
bash /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/tests/test_mock_pipeline.sh
bash /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/tests/test_region_extraction.sh
bash /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/tests/test_stranded_coverage.sh
```

Orkestreringstests bruger simulerede værktøjsresultater. Regionstesten
sammenligner regionsudtrækket med den tidligere metode på en syntetisk BAM
med rigtige samtools og sambamba; ingen af testene dokumenterer diagnostisk performance.

## Valgfri lokal SAI-10k-prædiktion

Tilføj `-sai-module /sti/til/sai10k_predictions.py` til analysekommandoen
for at aktivere den lokale Broad SpliceAI Lookup-implementering. Der kræves
ikke en præberegnet variantdatabase, og ingen variantoplysninger sendes til
en tjeneste. Modulet installeres separat fra Broad SpliceAI Lookup; den
afprøvede kildeversion er commit
`6412f6fdf5ae936bf732d88a7e3c9e5eeb9a272c` (MIT-licens). Kontrollér
modul og licens lokalt; et import-check alene beviser ikke, at filen
indeholder prædiktionsfunktionen:

```bash
MODULE=/drift/RNA/TotalRNA_Splice/Tools/SpliceAI-lookup/google_cloud_run_services/sai10k_predictions.py
wc -c "$MODULE"
conda run --no-capture-output -n KGA_Splice \
  python /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/scripts/sai_lookup.py \
  --module "$MODULE" --check
conda run --no-capture-output -n KGA_Splice \
  python /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/tests/test_sai_lookup.py "$MODULE"
```

Sidste test bruger kun syntetiske data og kontrollerer modulets resultat
for en deletion og en SNP. Standardmiljøet er `KGA_Splice`; angiv
`-sai-env MILJØNAVN` for et andet Conda-miljø med kompatibel Python.
Modulet, dets afhængigheder og `spliceai_env` skal være tilgængelige i
det miljø, hvor pipelinen køres.

SpliceAI beregnes en ekstra gang for det valgte GTF-transkript med
`-D 4999 -M 0`; det kan bruge mere tid og hukommelse end standardkørslen.
Transkriptet skal have gyldige exon- og CDS-annotationer. Simple indels
og SNP'er understøttes; komplekse substitutioner (både REF og ALT mere
end én base) angives som ikke understøttet. Ingen score eller variant
uden for transkriptet er **ikke** en negativ prædiktion.

Resultaterne er eksplorative, ikke observerede RNA-effekter. Broad-koden
får kun standard SpliceAI's fire delta-scores og positioner, ikke websidens
fulde positions-/REF-/ALT-scores, og anvender den lokale GTF. Resultatet
kan derfor afvige fra spliceailookup.broadinstitute.org. Frameshift
beregnes ud fra struktur, men stopkodon og proteinsekvens vurderes ikke.
Efterprøv prædiktioner i RNA-data, IGV og lokal validering. `STATUS.txt`
og `Report/index.html` angiver status, og `SAI10k/` gemmer input,
scores, `predictions.tsv`, detaljeret `predictions.json` og provenance
med modulets MD5. Ved utilgængelige scores kan nogle outputfiler mangle.

Kør derefter et testpar med kendt variant og kendt transkript. Opret først
en **tom, ny** output-overmappe. Erstat alle pladsholdere med korrekte
værdier; brug et kendt korrekt REF/ALT for GRCh38:

```bash
mkdir -p /drift/RNA/TotalRNA_Splice/TestRuns/ClinV1_test

Rscript /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1/scripts/ClinV1.R \
  -variant "CHR-POS-REF-ALT" \
  -gene GENE \
  -tx REFSEQ_TRANSCRIPT_VERSION \
  -p /drift/RNA/TotalRNA_Splice/TestRuns/ClinV1_test \
  -bam /drift/RNA/Drift \
  -s ID_NT,ID_PURO \
  -gt het \
  -refs /drift/RNA/TotalRNA_Splice/RefSample/TotalRefs
```

`-gt` er valgfri; `het`, `hmz`, `hemi` angiver patientens kendte
genomiske zygositet for varianten. Værdien verificeres ikke automatisk.
Uden `-gt` står ASE som ikke vurderet. Ved `het` og en biallelisk SNP
køres GATK ASEReadCounter for begge patientprøver, og REF-/ALT-tællinger
med rå CSV gemmes. Under 10 informative REF+ALT-reads mærkes resultatet
`low_depth_exploratory`; dette er en advarsel, ikke en valideret klinisk
grænse. Nul reads og manglende resultatrække vises særskilt. For
`hmz`/`hemi` eller en indel udføres ikke denne alleltest.

Der måles **kun den indtastede variant**, ikke andre heterozygote SNP'er
i genet. En SNP i samme gen kan give yderligere information, men uden
genomisk genotypeliste og fasning kan dens alleler ikke sikkert knyttes
til VUS. Alleltælling kan være nyttig for MAE/ASE, men er påvirket af
dækning, mappingsbias, duplikater, splejsning og NMD og viser ikke alene,
at en variant er årsag. Ingen binomial p-værdi eller automatisk
MAE-klassifikation anvendes uden lokal validering.

`-refs` er valgfri; udelad begge argumenter, hvis der ikke er relevante
reference-BAM'er. Alle BAM-filer i den valgte mappe inkluderes. En pooled
BAM er **én pool**, ikke individuelle biologiske kontroller, og blandede
behandlinger må ikke betragtes som en matchet kontrol. nt/puro-id'erne
skal være forskellige og kan indeholde bogstaver, tal, punktum, underscore
og bindestreg (første tegn skal være bogstav eller tal). Hvert ID skal matche
præcis én `ID.bam` eller `ID_*.bam`. For filerne `ID-B_*.bam` og
`ID-A_*.bam`, hvor B er nt og A er puro, anvendes `-s ID-B,ID-A`.
Behandlingsstatus skal bekræftes for hvert datasæt. Resultaterne ligger i én mappe:
`GENE_ID_NT_nt_ID_PURO_puro` under `-p`. Et eksisterende projekt
overskrives ikke. Ved fejl oprettes `STATUS.txt` med `INCOMPLETE`; brug en
ny output-overmappe ved næste forsøg.

## Output og fortolkning

* `STATUS.txt`: fuldførelsesstatus og begrænsninger. At kørslen er fuldført
  betyder ikke, at varianten kan vurderes.
* `GENE.vcf`, `GENE.ann.vep.vcf`, `GENE.spliceai.vcf`: variant og fælles
  lokale annotationer. VCF indeholder ikke en påstået patientgenotype.
* `variant_depth.tsv`: dækning for hver REF-base fra `samtools depth` med
  værktøjets standardfiltre. Den er ikke allelspecifik, og god dækning af
  positionen er ikke nok til at udelukke en splice-effekt.
* `Coverage/*.plus.bedgraph` og `*.minus.bedgraph`: dækning for reads med
  henholdsvis `XS:+` og fremadrettet alignment-flag og `XS:-` og
  omvendt alignment-flag, som i v6.4; minus-værdier er negative. Reads
  uden gyldigt XS-tag eller med modstridende flag udelades. Dette er ikke
  nødvendigvis transkriptstrand for et parret, strandet bibliotek.
  `samples.tsv` opgør antal inkluderede og udeladte alignments. Tomme
  bedgraphs udelades fra IGV-sessionen; dette betyder ikke nul ekspression.
* `ASE/allele_counts.tsv` og, hvor kørt, `*.ase_counts.csv`: alleltællinger
  for målvarianten og eksplicit status. GATK's filtrering adskiller sig fra
  `samtools depth`, så tællingerne kan ikke forventes at være identiske.
* `Report/index.html`: fælles HTML-rapport med `sashimi.pdf` fra sashimipy
  og to junction-søjlefigurer (rå og dybdenormaliseret); filerne skal
  forblive i samme mappe. Sashimipy
  bruger v6.4's fr-firststrand- og XS-antagelser og en grænse på 20 reads
  per junction; bibliotekets strandprotokol skal kontrolleres. Manglende
  buer er ikke bevis for fravær. Figuren er ikke en signifikansanalyse.
  SpliceAI DS-farver bruger de vejledende 0,2 / 0,5 / 0,8-grænser fra
  [SpliceAI](https://github.com/Illumina/SpliceAI); farve og tekst er ikke
  klinisk klassifikation. SAI-10k-effekter er kategoriske og har ikke en
  selvstændig 0-1-score; transkriptspecifikke SpliceAI-input vises separat,
  når de er tilgængelige. Deres farver er ikke validerede SAI-tærskler.
* `Junctions/junction_counts.tsv`: ét sæt observationer pr. eksakt genomisk
  intron `[intron_start_0, intron_end_0)` og strand, med tællinger pr.
  prøve; `[150,200)` svarer til intronbaser 151-200 i 1-baserede koordinater.
  `selected_transcript`, `other_gene_transcript`, `not_in_gene_annotation`,
  `strand_unknown` og `other_strand_or_contig` beskriver alene match til
  den anvendte GTF. Uannoteret er ikke automatisk patogent.
  `left_boundary_exon` og `right_boundary_exon` viser hhv. nærmeste exon-slut
  og exon-start for det valgte transkript ved intronets genomisk venstre og
  højre grænse. `exact` er et præcist match; `nearest` er et positionsmatch
  med absolut afstand i bp, ikke et bekræftet splice-site. Exoner nummereres
  efter GTF-koordinater i transkriptets 5'-til-3'-retning, ikke nødvendigvis
  efter GTF-attributten `exon_number`. Ved lige afstand vises alle kandidater.
  Også junctions med ukendt eller modsat strand sammenlignes positionsmæssigt;
  exonkolonnerne bekræfter ikke strand eller transkripttilhørsforhold. Andre
  transkripter kan have andre exonnumre og grænser. Kolonnerne findes også i
  rapportens top-junctiontabel.
  `junction_id` (J1, J2, ...) forbinder figur, tabel og TSV inden for én
  kørsel; ID'er kan ændre sig mellem kørsler. `variant_distance_bp` er den
  korteste afstand fra variantens 1-baserede REF-interval til første eller
  sidste intronbase; for insertion bruges VCF-ankerbasen. `near_variant`
  markerer afstand højst 100 bp på samme kromosom, også for ukendt/modsat
  strand. Andet kromosom giver `NA` i afstand. Grænsen er kun en visuel
  prioritering, ikke en valideret biologisk eller klinisk tærskel; relevante
  junctions kan ligge længere væk. Søjlefiguren viser top 20 samt op til 20
  ekstra nærliggende junctions efter afstand og dernæst reads, markeret med
  stjerne og rød tekst. Tabellen viser top 30 og alle nærliggende junctions.
  Kolonner med suffikset `_per_million_mapped` er rå junctiontællinger
  ganget med 1.000.000 og divideret med summen af mappede alignments i hele
  input-BAM ifølge `samtools idxstats`. Nævnerne findes i `samples.tsv` og
  de komplette idxstats-resultater i `Logs/`. Ved nul mappede alignments
  er værdierne `NA`, og prøven udelades fra den normaliserede figur. Dette
  er ikke PSI eller en statistisk test; sekundære/supplementære alignments,
  forskelle i filtre, genekspression, librarytype og poolens sammensætning
  kan påvirke sammenligningen. De to figurer bruger samme udvalg af junctions.
* `Junctions/*.bed`: rå regtools BED12; `BAM/`: indeksede udsnit;
  `samples.tsv`: input- og outputkobling; `Logs/`: kommando- og
  værktøjsoutput; `IGV/paired_IGV.xml`: fælles IGV-session. Region- og
  rå-BAM alignment-spor er som standard skjult, men deres ressourcer bevares;
  de kan vises igen i IGV. Coverage- og junctionspor skjules ikke.

Regtools køres med `-s XS -a 8 -m 50 -M 500000`. Manglende/ukorrekte
XS-tags kan gøre strandangivelsen ukendt eller forkert; verificér den mod
library-preparation og aligner før fortolkning. `.` , `*` og `?` bevares
som forskellige rå strandværdier, men klassificeres alle som ukendt strand
og matches ikke med `+` eller `-`. Junctions uden tilstrækkelig
anchor, uden for størrelsesgrænserne eller med utilstrækkelig dækning kan
mangle. Nul i tabellen betyder ingen junction registreret under disse
indstillinger, ikke ingen biologisk hændelse. Read counts er hverken PSI,
normaliseret ekspression eller statistiske test; ingen replikater er antaget.

XML indeholder containerens absolutte stier. IGV på en anden vært skal have
adgang til både regions-BAM'er og de originale patient- og kontrol-BAM'er
med indeks på tilsvarende stier, eller stierne skal oversættes. Patientspor
er blå i to nuancer for nt og puro; kontrolspor er grønne og har ukendt
behandling, medmindre den er verificeret separat. Variantsporet er rødt.
IGV-session, installation og biologisk resultat skal afprøves i brugerens
Linux-miljø. De syntetiske tests dokumenterer ikke diagnostisk performance.

Se [kohorteplanen](docs/KOHORTEPLAN.md) for OUTRIDER, FRASER, antal og
krav til validering.

Separate scripts til at oprette nt-/puro-BAM-pools og opbygge, udvide
og sammenligne en individuel junction-referencebank er beskrevet i
[REFERENCEBANK.md](docs/REFERENCEBANK.md). De er ikke automatisk koblet
til ClinV1 og er ikke klinisk validerede.
