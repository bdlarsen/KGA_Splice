# Docker: byg og kør ClinV1

ClinV1 er en **forsknings-/valideringskandidat**. Et succesfuldt build eller
en gennemført analyse er ikke klinisk validering eller godkendelse.
Docker-løsningen bygger oven på et eksisterende lokalt image og er ikke et
selvstændigt image fra bunden.

## Eksisterende snapshot: fælles kandidat til v6.4 og ClinV1

`Dockerfile.snapshot` er et alternativ til buildet nedenfor, når den aktuelle
container først er gemt med `docker commit`. Den bevarer snapshot-miljøerne
og gør `clinv1_regtools` tilgængelig på PATH uden at geninstallere eller
fjerne pakker. Buildet kontrollerer udvalgte afhængigheder, men beviser ikke
at alle moduler eller en fuld analyse fungerer. Snapshotten indeholder ikke
bind mounts; den skal bevares separat og kan indeholde følsomme lokale data.
Dette er ikke en fuldt reproducerbar opskrift på miljøet fra bunden.

På Docker-værten, med et eksisterende, testet snapshot:

```bash
BASE=/drift/RNA/TotalRNA_Splice
SNAPSHOT=kga-splice:container-backup-20261124-093511
docker image inspect "$SNAPSHOT" --format '{{.Id}}'
docker build --network none \
  --build-arg "BASE_IMAGE=$SNAPSHOT" \
  -f "$BASE/Pipeline/ClinV1/Dockerfile.snapshot" \
  -t kga-splice:combined-candidate \
  "$BASE/Pipeline/ClinV1"
```

Angiv det faktiske snapshot-tag og kontrollér image-ID før build. Kandidaten
ændrer ikke snapshotten, den kørende container eller `kga-splice:v6.3`.
Den indeholder ikke nye udgaver af scripts eller SAI-modulet: brug eksplicit
bind mounts for scripts, referencer, BAM og SAI, og bemærk at et projekt-mount
overskygger eventuelle scripts i image på samme sti. Test kandidaten og begge
pipelines med de tilsigtede mounts, før den eventuelt får et eksisterende tag.
Efter ændring af et tag fortsætter allerede oprettede containere med deres
oprindelige image; en ny container kræves. Bevar backup og foretag ikke
oprydning før validering. Ved runtime skal `regtools` stadig kontrolleres på
PATH efter eventuel Conda-aktivering.

## 1. Forudsætninger

På Linux-værten skal Docker virke, og projektet skal være tilgængeligt.
Eksemplet bruger det eksisterende image `kga-splice:v6.3`; vælg det image,
der faktisk indeholder dine fungerende miljøer. Kontrollér:

```bash
docker image inspect kga-splice:v6.3 --format '{{.Id}}'
```

Base-imaget skal have Conda under `/opt/conda`, miljøet `KGA_Splice` med
R, `xml2`, samtools, sambamba, bedtools, GATK og `sashimipy` fra `sashimi.py`,
samt fungerende `vep_env` og
`spliceai_env`. Dockerfilen tilføjer `regtools=1.0.0` i et separat
Conda-miljø, så den manglende `regtools` fra den tidligere kørsel ikke
kræver ændring af `KGA_Splice`.

Build kræver netværk til Conda-kanalerne, diskplads og en Docker-arkitektur,
som både base-imaget og pakkerne understøtter. Base-imaget, versioner i det
og Conda-afhængigheder er ikke låst fuldstændigt: et tag kan ændres, og
fremtidige builds kan få andre afhængigheder. Gem base-image-ID, image-ID
og relevante versionsoplysninger sammen med valideringsdokumentationen;
klinisk reproducerbarhed kræver yderligere fastlåsning og validering.

## 2. Byg image

Kør på **Linux-værten**, ikke inde i den gamle container:

```bash
cd /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1

docker build \
  --build-arg BASE_IMAGE=kga-splice:v6.3 \
  -t kga-splice:clinv1 \
  .
```

Er base-imaget fx `kga-splice:v6.4`, erstat `BASE_IMAGE` tilsvarende.
Hvis base-imaget kun findes i en eksisterende containers ændringer og ikke
som image, indeholder buildet **ikke** disse ændringer.

Dockerfilen tester værktøjer, Conda-miljøer og syntetiske tests under
build. Hvis en test eller afhængighed fejler, fejler buildet; kontrollér
den første fejl i outputtet. Særligt skal `vep_env` og `spliceai_env` findes
i base-imaget. Der downloades ikke FASTA, GTF, VEP-cache eller patientdata.
Disse filer skal leveres separat, og VEP-/SpliceAI-ressourcer skal passe
til GRCh38 og den installerede softwareversion.

```bash
docker image inspect kga-splice:clinv1 --format '{{.Id}}'
docker run --rm --network none kga-splice:clinv1 \
  bash -c 'command -v regtools && regtools --help'
```

Imaget indeholder et øjebliksbillede af `scripts/` og `tests/` under
`/opt/clinv1`. Ændringer på værten træder først i kraft efter et nyt
`docker build`. Ingen bred projektmappe mountes oven på disse scripts.

## 3. Start en container med begrænsede mounts

Eksemplet nedenfor bruger brugerens eksisterende reference-layout.
Erstat `BAM_DIR` med den **faktiske Linux-værtsmappe**, som indeholder
patientens BAM- og indeksfiler. Opret også en ny outputmappe på værten.
Alle stier til `--mount` skal eksistere, før containeren startes:

```bash
REF_ROOT=/drift/RNA/TotalRNA_Splice
BAM_DIR=/drift/RNA/REPLACE_WITH_BAM_DIRECTORY
OUT_DIR="$REF_ROOT/TestRuns/ClinV1_docker"
mkdir -p "$OUT_DIR"

docker run --rm -it --network none \
  --name clinv1-dev \
  --mount "type=bind,src=$BAM_DIR,dst=/data/bams,readonly" \
  --mount "type=bind,src=$REF_ROOT/Refs/references,dst=/refs/genome,readonly" \
  --mount "type=bind,src=$REF_ROOT/GRCh38_Reference.gtf,dst=/refs/GRCh38_Reference.gtf,readonly" \
  --mount "type=bind,src=$REF_ROOT/VEP_Refs,dst=/refs/vep,readonly" \
  --mount "type=bind,src=$OUT_DIR,dst=/data/output" \
  kga-splice:clinv1
```

`--network none` er valgt, fordi analysen og VEP skal køre lokalt/offline.
BAM, FASTA/GTF og cache er read-only, og kun outputmappen er skrivbar via
disse mounts. Containeren kan stadig skrive til sit eget filsystem. Docker
kører normalt som root her: kontrollér lokale adgangs- og ejerskabsregler,
og ret om nødvendigt ejerskab på output efter institutionens procedurer.
En read-only mount er ikke i sig selv en komplet sikkerheds- eller
adgangskontrolløsning.

Når du forlader shellen med `exit`, fjernes containeren pga. `--rm`.
Resultaterne i den mountede outputmappe bliver liggende. Andre ændringer
inde i containeren bevares ikke. Vil du medtage en pooled reference-BAM,
skal dens mappe mountes separat read-only, fx:

```text
--mount type=bind,src=/drift/RNA/TotalRNA_Splice/RefSample/TotalRefs,dst=/data/refs,readonly
```

Tilføj dette som en ekstra linje **før image-navnet** i `docker run`.
En pool er ikke uafhængige kontroller; `-refs` anvender samme pool til
begge konditioner, ikke én behandlingsmatchet pool til hver.

Ved valgfri SAI-10k skal det lokalt installerede Broad-modul og licensen
også mountes read-only. Tilføj før image-navnet:

```text
--mount type=bind,src=/drift/RNA/TotalRNA_Splice/Tools/SpliceAI-lookup/google_cloud_run_services,dst=/tools/spliceai-lookup,readonly
```

Tilføj derefter `-sai-module /tools/spliceai-lookup/sai10k_predictions.py`
til analysekommandoen. Modulet er ikke inkluderet i imaget. Efter ændringer
af ClinV1-scripts skal imaget bygges igen. Krav, syntetisk test af modulet
og begrænsninger er beskrevet i [README](../README.md).

## 4. Test inde i containeren

```bash
cd /opt/clinv1
Rscript tests/test_core.R
Rscript tests/test_report.R
Rscript tests/test_sai.R
Rscript tests/test_report_pipeline.R
bash tests/test_mock_pipeline.sh
bash tests/test_region_extraction.sh
bash tests/test_stranded_coverage.sh
Rscript tests/test_reference_pool.R
Rscript tests/test_reference_bank.R
bash tests/test_reference_workflows.sh
```

Tests er syntetiske og beviser ikke biologisk eller klinisk performance.
`test_mock_pipeline.sh` fremkalder en forventet værktøjsfejl og skal slutte
med `Mock pipeline tests passed`. Reference-workflowtesten fremkalder også
en forventet afvisning af patienten som sin egen kontrol og skal slutte med
`Reference workflow integration tests passed`.

## 5. Kør en analyse

Alle følgende værdier er pladsholdere. Brug en kendt GRCh38-variant og
angiv først ID uden puromycin, derefter ID med puromycin:

```bash
Rscript /opt/clinv1/scripts/ClinV1.R \
  -variant "CHR-POS-REF-ALT" \
  -gene GENE \
  -tx REFSEQ_TRANSCRIPT_VERSION \
  -p /data/output \
  -bam /data/bams \
  -s ID_NT,ID_PURO \
  -gt het \
  -g /refs/genome/genome.fa \
  -gtf /refs/GRCh38_Reference.gtf \
  -vep /refs/vep
```

Outputprojektet må ikke eksistere i forvejen. Fejl efter projektoprettelse
giver en `INCOMPLETE`-status; kontrollér konsol og `Logs/`, og anvend en ny
outputmappe ved næste kørsel. En fejl i preflight kan forekomme, før
projektmappen oprettes. Se [README](../README.md) for output og
fortolkningsbegrænsninger.

## 6. Referencebank i container

Pool- og bank-scripts ligger også i `/opt/clinv1/scripts/`. Beskrivelserne
og manifestformaterne findes i [REFERENCEBANK.md](REFERENCEBANK.md). Brug
de **stier, der findes inde i containeren**, i BAM-lister og manifest.
Mount individuelle BAM-mapper read-only og en separat bank-/poolmappe
read-write. En bank gemmer inputstier til kontrol ved senere tilføjelse;
brug derfor konsekvent samme containerstier på tværs af kørsler.

## 7. Hvis noget fejler

* `regtools: not found`: kontrollér, at den nye image-tag bruges, og at
  shellen ikke efterfølgende har ændret PATH; `command -v regtools`.
* `vep_env` eller `spliceai_env` mangler: vælg et base-image med de
  fungerende miljøer, eller opbyg disse i et separat godkendt base-image.
* `mount source path does not exist`: kontrollér stier på **Docker-værten**;
  `--mount` opretter ikke manglende kilder.
* Manglende VEP-/SpliceAI-annotation er ikke det samme som en negativ
  prædiktion. Gennemgå `STATUS.txt` og logs, og kontrollér cache/reference.
* IGV-sessionen indeholder containerstier. IGV på en anden vært skal have
  adgang til filerne med tilsvarende stier, eller stierne skal oversættes.
