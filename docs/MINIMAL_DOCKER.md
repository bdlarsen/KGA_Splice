# ClinV1: minimalt Docker-kandidatimage og testplan

Dette image er et selvstændigt **forsknings-/valideringskandidatimage** til
ClinV1, ikke til v6.4. Det erstatter eller ændrer ikke eksisterende images.
Der installeres R/xml2, samtools, sambamba, bedtools, regtools, GATK,
sashimipy og separate miljøer til VEP og SpliceAI. Broad SAI-10k-modulet og
referencer mountes separat. Det indeholder ikke IGV Desktop/Xvfb eller HGVS:
IGV-snapshot er fortsat eksperimentelt og HGVS-input ikke implementeret.

Dette er en opskrift fra bunden, men ikke et fuldt versionslåst eller
valideret build. Base-image-tag og indirekte afhængigheder kan ændre sig.
VEP-cache, FASTA, GTF, reference-BAM, modeldata og softwareversioner skal
passe sammen og valideres lokalt. VEP er specificeret til release 115;
cache skal være en kompatibel RefSeq-cache for GRCh38. SpliceAI bruger sin
indbyggede grch38-annotation. Dokumentér de faktiske versioner og image-ID.
Build og package-download kræver netværk; patientanalyse køres uden netværk.
Vær opmærksom på SpliceAI- og øvrige pakkers licenser og lokale brugsregler.
Base-miljøets inkompatible telemetry-plugin fjernes målrettet under build:
Conda-pakkemetadata kontrolleres først for afhængigheder, og buildet stopper,
hvis pluginet mangler eller en base-pakke afhænger af det. Tvungen fjernelse
omgår Condas solver og er ikke en generel garanti for kompatibilitet.

## 1. Byg på Docker-værten

Vent til eventuelle igangværende kørsler er afsluttet. Kontrollér plads i
Dockers dataområde; bygning kan kræve betydeligt mere disk end slutimaget.
Der er kun rapporteret cirka 12 GB ledigt i containerens overlay, så sørg
for tilstrækkelig plads før build. Opdater ikke eksisterende image-tags.

```bash
BASE=/drift/RNA/TotalRNA_Splice
CTX="$BASE/Pipeline/ClinV1"
IMAGE=kga-splice:clinv1-minimal-candidate

docker info --format 'Docker data: {{.DockerRootDir}}'
df -h /var/lib/docker "$BASE"
docker build --progress=plain -f "$CTX/Dockerfile.minimal" -t "$IMAGE" "$CTX"
```

Hvis Docker-data ikke ligger under `/var/lib/docker`, kontrollér pladsen på
den faktiske sti. `latest` og de ulåste afhængigheder kræver netværk og kan
give forskellige resultater senere; behold og dokumentér det testede image.
Stop ved første buildfejl; ret ikke det gamle image for at få buildet igennem.

### Alternativ: Docker Compose

Compose-filen bygger samme kandidatimage. Kør på værten; build kræver netværk,
men containere startet af denne service har `network_mode: none`.
Et fast værtsnavn og en lokal loopback-post sikrer lokalt navneopslag for GATK
uden at åbne for netværk.

```bash
BASE=/drift/RNA/TotalRNA_Splice
CTX="$BASE/Pipeline/ClinV1"
COMPOSE="$CTX/compose.minimal.yaml"
IMAGE=kga-splice:clinv1-minimal-candidate

docker compose version
docker compose -f "$COMPOSE" config
docker compose --progress=plain -f "$COMPOSE" build
docker image inspect "$IMAGE" --format 'ID={{.Id}} Size={{.Size}} bytes'
docker compose -f "$COMPOSE" run --rm -T --no-deps clinv1 \
  bash /opt/clinv1/tests/test_minimal_image.sh
```

Start ikke patientanalyse før de syntetiske tests og de efterfølgende
valideringskontroller er gennemført. `run --rm` fjerner kun den nye
testcontainer; image og eksisterende containere berøres ikke.

Til analyse kan Compose også bruges med eksplicitte mounts. Definér de
værtsstier, der vises i afsnit 4, og start derefter:

```bash
if [[ -d "$BAM_DIR" && -d "$BASE/Refs/references" \
   && -f "$BASE/GRCh38_Reference.gtf" && -d "$BASE/VEP_Refs" \
   && -d "$BASE/RefSample/TotalRefs" && -d "$SAI" && -d "$OUT_DIR" ]]; then
  docker compose -f "$COMPOSE" run --rm --no-deps \
    -v "$BAM_DIR:/data/bams:ro" \
    -v "$BASE/Refs/references:/refs/genome:ro" \
    -v "$BASE/GRCh38_Reference.gtf:/refs/reference.gtf:ro" \
    -v "$BASE/VEP_Refs:/refs/vep:ro" \
    -v "$BASE/RefSample/TotalRefs:/data/refs:ro" \
    -v "$SAI:/tools/sai:ro" \
    -v "$OUT_DIR:/data/output" \
    clinv1 bash
else
  echo "FEJL: En eller flere værtsstier mangler. Containeren blev ikke startet."
fi
```

Kontrollér først at alle værtsstier findes; `-v` kan ellers oprette en
mappe på værten. Følg herefter analysekommandoen og testplanen i afsnit 4
og 5. Mounts er kun tilføjet ved det konkrete `run`, ikke i Compose-filen.

## 2. Kontroller imaget uden patientdata

På værten:

```bash
docker image inspect "$IMAGE" --format 'ID={{.Id}} Size={{.Size}} bytes'
docker run --rm --network none --hostname clinv1-local \
  --add-host clinv1-local:127.0.0.1 "$IMAGE" \
  bash /opt/clinv1/tests/test_minimal_image.sh
```

Forvent `Minimal image synthetic tests passed.` til sidst og exitkode 0.
Testen omfatter simulerede pipelinekørsler, parser-, rapport-, VCF-, enkelt-
og parprøvetests, samt syntetiske tests af reelle samtools, sambamba,
bedtools og regtools. Den kontrollerer ikke en fuld rigtig VEP-/SpliceAI-
analyse, IGV-rendering eller biologisk performance. `IGV snapshot failed`
i testoutput er forventet, når fejlhåndtering testes.

Gem desuden versionslister sammen med jeres valideringsdokumentation:

```bash
docker run --rm --network none --hostname clinv1-local \
  --add-host clinv1-local:127.0.0.1 "$IMAGE" bash -c \
  'conda list -n KGA_Splice && conda list -n vep_env && conda list -n spliceai_env'
```

## 3. Kontrollér SAI-modul med syntetiske data

På værten; stien skal eksistere og være den kontrollerede Broad-version:

```bash
SAI="$BASE/Tools/SpliceAI-lookup/google_cloud_run_services"
test -s "$SAI/sai10k_predictions.py" || echo "FEJL: SAI-modul mangler"
docker run --rm --network none --hostname clinv1-local \
  --add-host clinv1-local:127.0.0.1 \
  --mount "type=bind,src=$SAI,dst=/tools/sai,readonly" "$IMAGE" \
  conda run --no-capture-output -n KGA_Splice \
  python /opt/clinv1/tests/test_sai_lookup.py /tools/sai/sai10k_predictions.py
```

Kontrollér både exitkode og testresultat. Ingen patientvarianter bruges her.
Dette verificerer ikke de rigtige lokale VEP-/SpliceAI-referencer.

## 4. Fuld test på kendt prøvepar

Vælg en allerede godkendt lokal valideringsprøve med kendt facit. Sæt
værdierne nedenfor på **værten**; `BAM_DIR` skal indeholde begge patient-BAM
og indeks direkte i samme mappe. Opret et nyt outputområde på værten; undgå
at genbruge en tidligere projektmappe. `-gt het` er kun et eksempel og må
kun anvendes, hvis genotypen er uafhængigt bekræftet.

```bash
BASE=/drift/RNA/TotalRNA_Splice
IMAGE=kga-splice:clinv1-minimal-candidate
BAM_DIR=/sti/til/valideringsbam
OUT_DIR="$BASE/TestRuns/ClinV1_minimal_validering"
SAI="$BASE/Tools/SpliceAI-lookup/google_cloud_run_services"
mkdir -p "$OUT_DIR"

docker run --rm -it --network none --hostname clinv1-local \
  --add-host clinv1-local:127.0.0.1 \
  --mount "type=bind,src=$BAM_DIR,dst=/data/bams,readonly" \
  --mount "type=bind,src=$BASE/Refs/references,dst=/refs/genome,readonly" \
  --mount "type=bind,src=$BASE/GRCh38_Reference.gtf,dst=/refs/reference.gtf,readonly" \
  --mount "type=bind,src=$BASE/VEP_Refs,dst=/refs/vep,readonly" \
  --mount "type=bind,src=$BASE/RefSample/TotalRefs,dst=/data/refs,readonly" \
  --mount "type=bind,src=$SAI,dst=/tools/sai,readonly" \
  --mount "type=bind,src=$OUT_DIR,dst=/data/output" \
  "$IMAGE" bash
```

Inde i containeren, erstat alle pladsholdere:

```bash
Rscript /opt/clinv1/scripts/ClinV1.R \
  -variant 'CHR-POS-REF-ALT' -gene GENE -tx NM_VERSION \
  -p /data/output -bam /data/bams -s ID_NT,ID_PURO -gt het \
  -g /refs/genome/genome.fa -gtf /refs/reference.gtf -vep /refs/vep \
  -refs /data/refs -sai-module /tools/sai/sai10k_predictions.py
```

For en enkelt patientprøve bruges `-s ID`; behandling sættes til ukendt.
Indekseret FASTA og GATK-dictionary skal ligge ved siden af FASTA. IGV XML
indeholder containerstier, der muligvis skal oversættes til IGV-værten.
Der må ikke tilføjes netværk for at få annotering til at fungere.

## 5. Acceptkriterier før anvendelse

* Gem buildlog, image-ID, pakkelister, versions- og referenceoplysninger.
* Kontrollér `STATUS.txt`, logs, REF-kontrol, VEP og SpliceAI-annotation,
  SAI-status, `samples.tsv`, dækning, junctiontællinger, ASE og rapport.
  Lav eller manglende dækning skal fremgå og må ikke fortolkes som negativt.
* Sammenlign samme input med tidligere gennemgået kandidatimage og kendt
  facit, både rå data og visning. Forskel i software/referenceversioner kan
  påvirke output; undersøg forskelle før accept. Kontrollér nt/puro-identitet.
* Gentag for enkeltprøve, heterozygot SNP, indel, nul/lav dækning og en
  forventet fejl, herunder REF-mismatch og manglende BAM-indeks.
  Kontroller at fejl stopper og ikke producerer succesformet status.
* Åbn IGV-sessionen manuelt og kontrollér BAM, coverage og junctiontracks.
  Den automatiske IGV-SVG-funktion er ikke valideret eller installeret her.
* Udvid med lokale positive og negative cases på begge strands, flere gener,
  transkriptversioner, bibliotekstyper og relevante referencekontroller.
  Afvigelser, begrænsninger, licenser og resultatgodkendelse skal håndteres
  efter lokal protokol. Et vellykket build er ikke klinisk validering.

Kun efter separat godkendelse kan kandidaten tages i brug eller få et andet
tag. Bevar det nuværende fungerende image og referencer; slet intet som del
af denne test. Outputs indeholder følsomme oplysninger og skal opbevares efter
lokale regler.
