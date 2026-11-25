# Git og GitHub: trin for trin

ClinV1 er klargjort som et lokalt Git-repository. Denne guide udfører
**ingen upload**. Sørg først for, at organisationen tillader, at kode og
dokumentation deles på den valgte GitHub-konto. Vælg et privat repository,
men vær opmærksom på, at privat GitHub stadig er en ekstern tjeneste.

**Upload aldrig** patient-ID'er, varianter, BAM/FASTQ, VCF, resultater,
referencebank, logs, lokale manifestfiler, credentials eller reference-cache.
Dokumentation, kode, billedfiler, commit-historik og filnavne skal også
gennemgås manuelt. `.gitignore` og `.dockerignore` mindsker risikoen, men
garanterer ikke anonymisering eller fravær af følsomme oplysninger.
Send heller ikke Docker-image eller base-image til en registry uden separat
kontrol og godkendelse; et image kan indeholde mere end dette repository.

## 1. Gå til projektet på Linux

```bash
cd /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1
git status
```

Hvis mappen er en kopi uden `.git`, og `git status` siger, at den ikke er
et repository, kan du efter kontrol af stien oprette ét:

```bash
git init -b main
```

Kør ikke `git init` igen, hvis repositoryet allerede findes. På et
netværksdrev kan Linux og Windows se forskellige Git-konfigurationer;
kontrollér resultatet fra det miljø, du bruger til upload.

Hvis Git på det delte drev viser `detected dubious ownership`, bør
ejerskab og adgang afklares med administratoren. Er den præcise mappe
bekræftet betroet, kan en enkelt Linux-kommando køres med en midlertidig
undtagelse, uden at ændre global Git-konfiguration:

```bash
git -c safe.directory="$(pwd -P)" status --short --untracked-files=all
```

Undtagelsen gælder kun den ene kommando. Det løser ikke nødvendigvis en
anden ejerskabs- eller adgangsfejl. Aktivér ikke blindt tillid til alle
repositories med `safe.directory=*`.

## 2. Gennemgå filerne før staging

```bash
git status --short --untracked-files=all
```

Kun kode, tests, Dockerfil og dokumentation skal med. Ignorerede filer kan
du se med `git status --short --ignored`. Ignorerede filer skal ikke
force-addes med `git add -f`. Gennemgå alle kilder og eksempler manuelt
for følsomme oplysninger — også oplysninger, der ikke har kendte filendelser.
Git-ignore-reglerne er en allowlist, ikke en datasikkerhedsscanner.
Store referencefiler og patientdata må ikke placeres i repositoryet, heller
ikke via Git LFS.

## 3. Stage de godkendte filer og kontrollér dem

```bash
git add .gitattributes .gitignore .dockerignore Dockerfile README.md scripts tests docs
git diff --cached --name-only
git diff --cached --check
git diff --cached
```

Kontrollér både filnavne og indhold. Hvis et følsomt dokument ved en fejl
er staged, fjern det fra staging med `git restore --staged -- STIEN` før
commit, og håndtér selve filen efter lokale regler. Fortsæt ikke, hvis der
er tvivl. Når følsomt indhold først er committet eller pushet, er fjernelse
fra seneste version ikke det samme som at fjerne det fra historikken.

## 4. Opret første commit

Git skal have et godkendt navn og en e-mailadresse. Undgå en privat eller
uønsket adresse i et delt repository; kontrollér organisationens politik.
Efter manuel godkendelse:

```bash
git commit -m "Initial ClinV1 pipeline candidate"
git status
git log -1 --oneline
```

Ved manglende Git-identitet skal du konfigurere den lokalt i dette
repository med `git config user.name` og `git config user.email` ifølge
organisationens vejledning, og derefter gentage commit.

## 5. Opret et tomt privat GitHub-repository

1. Log på den GitHub-konto/organisation, der er godkendt af din arbejdsplads.
2. Vælg **New repository**, angiv et passende navn og vælg **Private**.
3. Opret det som et **tomt** repository uden automatisk README, `.gitignore`
   eller licens, så historikkerne ikke divergerer.
4. Vælg selv passende licens efter aftale med rettighedshaver; denne pakke
   tilføjer ikke en licens. Uden eksplicit licens følger ingen generel
   tilladelse til viderebrug.
5. Kopiér repositoryets HTTPS- eller SSH-adresse fra GitHub.

## 6. Tilknyt og upload — kun efter godkendelse

Erstat URL'en med din egen; eksemplet er en pladsholder:

```bash
git remote add origin https://github.com/OWNER/REPOSITORY.git
git remote -v
git branch --show-current
git push -u origin main
```

Brug det faktiske branchenavn, hvis det ikke er `main`. GitHub kræver en
godkendt autentificeringsmetode for HTTPS/SSH. Skriv ikke tokens i URL,
kommandohistorik eller filer. Hvis `origin` allerede findes, undersøg først
`git remote -v`; overskriv ikke en eksisterende remote uden at forstå den.

## 7. Senere ændringer

```bash
git status --short --untracked-files=all
# Gennemgå ændringer, filnavne og nye filer manuelt.
git add .gitattributes .gitignore .dockerignore Dockerfile README.md scripts tests docs
git diff --cached --name-only
git diff --cached --check
git diff --cached
git commit -m "Describe the approved change"
git push
```

Pakkede Docker-images, BAM-filer, referencebank og resultater holdes
uden for Git. Ved utilsigtet offentliggørelse: stop upload/deling og følg
institutionens procedurer for informationssikkerhed. Overvej separat
gennemgang af hele commit-historikken før deling.
