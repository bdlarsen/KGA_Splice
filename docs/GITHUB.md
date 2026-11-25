# Git og GitHub: backup og daglig brug

Denne guide bruges på Linux-værten til at vedligeholde den oprindelige
ClinV1-mappe og dens private GitHub-backup. Kommandoerne kører ikke automatisk:
gennemgå altid indholdet, før du selv opretter et commit eller uploader.

Git gemmer versioner **lokalt**; GitHub er den eksterne kopi. `git add`
vælger ændringer til næste version, `git commit` opretter versionen lokalt,
og `git push` uploader commits. Et privat repository er stadig en ekstern
tjeneste: kontrollér, at anvendelsen er godkendt efter lokale regler.

## Hvad skal og skal ikke sikkerhedskopieres?

Repositoryet er til kode, Dockerfiler, tests og dokumentation. BAM/FASTQ,
VCF, referencebank, reference-cache, resultater, logs, lokale manifestfiler
og adgangsoplysninger skal ikke uploades. Prøve-ID'er er i dette projekt
oplyst at være anonymiserede; vurder alligevel indhold og eksempler efter
lokale regler. Fjern eventuelle personnavne fra scripts, før de deles.

Gennemgå også filnavne, dokumentation, kode og commit-historik manuelt.
`.gitignore` er en allowlist, ikke en sikkerhedsscanner; brug ikke
`git add -f` til at omgå den. Et commit kan gemme data i historikken, selv
om filen senere slettes. Stop og følg lokale sikkerhedsprocedurer, hvis
følsomt indhold allerede er committet eller uploadet.

Docker-images og referencedata er ikke del af denne kodebackup. En klon
er derfor ikke i sig selv en komplet, reproducerbar driftsinstallation.
Upload ikke images til en registry uden separat kontrol og godkendelse.

## 1. Start på hosten

Kør dette i hver ny terminalsession:

```bash
cd /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1
SAFE="$(pwd -P)"

git -c safe.directory="$SAFE" status --short --untracked-files=all
git -c safe.directory="$SAFE" branch --show-current
git -c safe.directory="$SAFE" remote -v
```

`SAFE` er en midlertidig undtagelse for netop denne mappe, hvis Git på det
delte drev viser `detected dubious ownership`. Brug den kun, når du har
bekræftet, at mappen er betroet. Undtagelsen gælder kun kommandoen, den
er angivet på; den ændrer ikke global konfiguration. Brug ikke
`safe.directory=*`. En anden adgangs- eller låsefejl skal undersøges separat.

Repositoryet er allerede oprettet: kør ikke `git init` igen. `origin` er
normalt navnet på GitHub-kopien. Kontrollér adressen og at GitHub viser
repositoryet som **Private**, før du uploader.

## 2. Se hvad der er ændret

Brug disse kommandoer efter arbejde på pipelinen:

```bash
git -c safe.directory="$SAFE" status --short --untracked-files=all
git -c safe.directory="$SAFE" --no-pager diff
git -c safe.directory="$SAFE" --no-pager diff --cached
```

`status` viser både nye, ændrede og slettede filer. Typiske markeringer er
`??` for en ny fil, `M` for en ændring og `D` for en sletning. Første
kolonne angiver staged ændringer, anden kolonne ændringer i arbejdskopien.
`diff` viser ikke-stagede ændringer i filer, som Git allerede følger;
nye `??`-filer skal åbnes og gennemgås særskilt. `diff --cached` viser
det, der vil indgå i næste commit.

Brug `git -c safe.directory="$SAFE" status --short --ignored` til at se
ignorerede filer. Tilføj ikke data bare fordi de kan ses i status.

## 3. Vælg ændringer til et commit

Når du kun vil gemme bestemte godkendte filer, er det sikrest at angive
dem eksplicit. Eksempel for denne guide og launcheren:

```bash
git -c safe.directory="$SAFE" add -- docs/GITHUB.md run_clinv1.sh
```

Hvis alle de relevante ændringer er gennemgået og skal med, kan de kendte
kode- og dokumentationsområder stages samlet:

```bash
git -c safe.directory="$SAFE" add -- \
  .gitattributes .gitignore .dockerignore \
  Dockerfile Dockerfile.minimal Dockerfile.snapshot \
  compose.minimal.yaml run_clinv1.sh README.md scripts tests docs
```

Staging tager et øjebliksbillede. Hvis du redigerer filen igen bagefter,
skal du køre `git add` på den igen for at få den nyeste udgave med.
Kontrollér også sletninger i status. Hvis en angivet fil er slettet, kan
du efter kontrol stage sletningen med eksempelvis
`git -c safe.directory="$SAFE" add -u -- docs/GITHUB.md`; den brede
kommando ovenfor kan fejle, hvis en eksplicit angivet rodfil mangler.

## 4. Kontrollér det staged indhold

```bash
git -c safe.directory="$SAFE" diff --cached --name-status
git -c safe.directory="$SAFE" diff --cached --check
git -c safe.directory="$SAFE" --no-pager diff --cached
```

Den første viser de valgte filnavne og ændringstyper, den anden kontrollerer
blandt andet whitespace-fejl, og den sidste viser selve indholdet. En tom
udskrift fra `--check` er ikke en sikkerhedsgodkendelse. Gennemgå indholdet
for personnavne, adgangsoplysninger og data, før du fortsætter. Del ikke
følsomt diff-indhold i chat eller andre eksterne tjenester.

Hvis en allerede versionsstyret fil er staged ved en fejl, kan den fjernes
fra staging uden at ændre arbejdskopien:

```bash
git -c safe.directory="$SAFE" restore --staged -- docs/GITHUB.md
```

Er filen ny og endnu ikke i seneste commit, kan dette i stedet bruges:

```bash
git -c safe.directory="$SAFE" rm --cached -- sti/til/ny-fil
```

Kontrollér stien først; `--cached` bevarer den lokale fil. Bekræft bagefter
med `git status` og `git diff --cached --name-status`.

## 5. Gem versionen lokalt: commit

Kontrollér først afsenderidentiteten; navn og e-mail gemmes i historikken
og deles med repositoryet:

```bash
git -c safe.directory="$SAFE" config user.name
git -c safe.directory="$SAFE" config user.email
```

Hvis Git mangler identitet, kan en godkendt identitet sættes **lokalt i
repositoryet** efter organisationens regler:

```bash
git -c safe.directory="$SAFE" config --local user.name "Godkendt navn"
git -c safe.directory="$SAFE" config --local user.email "godkendt@example.org"
```

Når staging er kontrolleret, opret et commit med en beskrivende tekst:

```bash
git -c safe.directory="$SAFE" commit -m "Update GitHub backup guide"
git -c safe.directory="$SAFE" status --short --untracked-files=all
git -c safe.directory="$SAFE" --no-pager log -5 --oneline
```

Et commit er stadig kun lokalt. Ændringer, der ikke var staged, er ikke
med. `log` viser de seneste commits; det er nyttigt til at kontrollere,
at versionen blev oprettet.

## 6. Upload godkendte commits

Når `origin` allerede findes, og den aktuelle branch følger GitHub-branchen:

```bash
git -c safe.directory="$SAFE" remote -v
git -c safe.directory="$SAFE" branch -vv
git -c safe.directory="$SAFE" push
```

Ved første push af en lokal `main`-branch uden upstream bruges i stedet:

```bash
git -c safe.directory="$SAFE" push -u origin main
```

`-u` opretter forbindelsen til fjernbranchen, så senere `git push` normalt
er nok. Brug det faktiske branchenavn, hvis det ikke er `main`. Kontrollér
derefter filerne på GitHub. Kør ikke `remote add origin` igen, når `origin`
allerede findes, og brug ikke force-push for at omgå fejl.

### Godkendelse med HTTPS

GitHub accepterer ikke kontoadgangskoden som adgangskode til Git over HTTPS.
Et fine-grained personal access token kan oprettes i GitHubs indstillinger
med adgang kun til det relevante repository og **Contents: Read and write**.
Metadata er normalt automatisk **Read-only**; Administration er ikke
nødvendig til almindeligt push. Vælg en passende udløbsdato.

Når Git spørger, indtast GitHub-brugernavn og indsæt tokenet ved
`Password`-prompten. Del aldrig tokenet, og skriv det ikke i en kommando,
URL, fil eller chat. En credential helper kan gemme det: kontrollér lokal
politik og værtsopsætning, særligt på en delt maskine. Alternativt kan
godkendt SSH eller en godkendt credential manager anvendes.

### Alternativ: SSH med en deploy key

En deploy key er knyttet til ét GitHub-repository. Med skriveadgang kan
den bruges til push uden et GitHub-token. En nøgle med passphrase kan
stadig kræve indtastning; en SSH-agent kan huske den i en session, hvis
det er tilladt efter lokale regler.

På en delt server skal brug af nøglen godkendes efter lokale regler.
Nøglen ligger i den aktuelle brugers hjemmemappe; bruges `root`, kan andre
med root-adgang potentielt også bruge den. Deploy keys har ikke en
automatisk udløbsdato som et tidsbegrænset token, så fjern adgangen på
GitHub, når den ikke længere er nødvendig.

**1. Opret nøglen på hosten**

```bash
ssh-keygen -t ed25519 -f "$HOME/.ssh/id_ed25519_clinv1" \
  -C "ClinV1 backup"
```

Hvis nøglen allerede findes, må den ikke overskrives. Brug den eksisterende
nøgle efter kontrol, eller vælg et andet filnavn og tilpas de efterfølgende
kommandoer. Vælg en passphrase, når du bliver spurgt.

**2. Tilføj den offentlige nøgle på GitHub**

```bash
cat "$HOME/.ssh/id_ed25519_clinv1.pub"
```

Kopiér hele linjen, der starter med `ssh-ed25519`, inklusive eventuel
kommentar til sidst. På GitHub vælges **KGA_Splice → Settings → Deploy
keys → Add deploy key**. Indsæt linjen i **Key**, angiv et navn, og vælg
**Allow write access**, før nøglen tilføjes. Uden skriveadgang kan den
ikke bruges til push.

Filen uden `.pub` er den private nøgle: den må aldrig uploades eller
deles. En offentlig nøgle skal være én linje i OpenSSH-format; kopier
hverken terminalprompten eller et linjeskift midt i nøglen.

**3. Test forbindelsen**

```bash
ssh -i "$HOME/.ssh/id_ed25519_clinv1" \
  -o IdentitiesOnly=yes -T git@github.com
```

Første gang kan SSH bede om godkendelse af GitHubs værtsnøgle. Kontrollér
fingeraftrykket mod [GitHubs officielle SSH-fingeraftryk](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/githubs-ssh-key-fingerprints)
før accept. En vellykket test bekræfter autentificering og siger, at GitHub
ikke giver shell-adgang; kommandoen afsluttes alligevel normalt med kode 1.
Ved anden fejl: stop og undersøg den, før du ændrer repositoryets remote.

**4. Brug SSH til dette repository**

Efter vellykket test kontrolleres først den eksisterende remote og eventuel
SSH-konfiguration:

```bash
cd /drift/RNA/TotalRNA_Splice/Pipeline/ClinV1
SAFE="$(pwd -P)"
git -c safe.directory="$SAFE" remote -v
git -c safe.directory="$SAFE" config --local --get core.sshCommand
```

Hvis `origin` er den tilsigtede GitHub-backup, og der ikke er en eksisterende
`core.sshCommand`, som skal bevares, kan SSH vælges lokalt for repositoryet:

```bash
git -c safe.directory="$SAFE" config --local core.sshCommand \
  "ssh -i $HOME/.ssh/id_ed25519_clinv1 -o IdentitiesOnly=yes"

git -c safe.directory="$SAFE" remote set-url origin \
  git@github.com:bdlarsen/KGA_Splice.git

git -c safe.directory="$SAFE" remote -v
git -c safe.directory="$SAFE" push
```

Undersøg eksisterende SSH-konfiguration først, hvis kommandoen viser en
værdi. `remote set-url` ændrer den eksisterende `origin`; kontrollér altid
adressen før push. Fortsæt ikke ved fejl. Skift ikke til en anden Linux-bruger
uden at kontrollere nøgle, hjemmemappe og adgang.

## 7. Hent ændringer fra GitHub

Brug kun dette, hvis en anden godkendt ændring faktisk skal ind på hosten.
At hente og indarbejde ændringer kan påvirke den pipeline, der køres i drift;
gennemgå og validér ændringerne efter lokale procedurer.

```bash
git -c safe.directory="$SAFE" status --short --untracked-files=all
git -c safe.directory="$SAFE" fetch origin
git -c safe.directory="$SAFE" --no-pager log --oneline --decorate --left-right main...origin/main
```

`fetch` henter information uden at ændre arbejdskopiens filer. I loggen
angiver `<` lokale commits og `>` commits på GitHub. Eksemplet forudsætter
branchen `main`; kontrollér `branch -vv`, hvis den hedder noget andet.
Når arbejdskopien er ren, ændringerne er gennemgået, og branchen har korrekt
upstream, kan en ren fast-forward udføres:

```bash
git -c safe.directory="$SAFE" pull --ff-only
```

Hvis kommandoen afvises, eller lokale og eksterne commits er divergeret,
stop og undersøg situationen; brug ikke `reset --hard` eller force-push.

## 8. Nyt, tomt repository (kun ved en ny backup)

Den eksisterende ClinV1-backup kræver ikke dette trin. Hvis en ny backup
skal oprettes, opret først et tomt **Private** repository på en godkendt
GitHub-konto. Vælg ikke automatisk README, `.gitignore` eller licens.
Vælg en eventuel licens efter aftale med rettighedshaver.

Kontrollér `git remote -v`. Kun hvis `origin` ikke allerede findes, kan
den tilføjes med den rigtige adresse uden adgangsoplysninger:

```bash
git -c safe.directory="$SAFE" remote add origin https://github.com/OWNER/REPOSITORY.git
git -c safe.directory="$SAFE" remote -v
```

Erstat pladsholderne, kontrollér adressen og brug første-push-kommandoen
ovenfor efter gennemgang og commit. Overskriv aldrig en eksisterende
remote uden først at forstå, hvad den bruges til.
