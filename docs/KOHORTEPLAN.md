# Implementeringsplan: kohortebaseret RNA-analyse

## Formål og afgrænsning

Målet er at supplere målrettet vurdering af VUS med en mere systematisk,
genomomfattende søgning efter ekspressions- og splice-outliers. OUTRIDER
vurderer genekspression; FRASER vurderer afvigende splice-junction-/intron-
brug. De er ikke direkte tests for, om netop den undersøgte variant er
årsag, og kan ikke erstatte vurdering af dækning, transkript, prøveidentitet,
genotype, relevante kontroller og orthogonal bekræftelse.

Én nt- og én puroprøve fra samme patient er forskellige betingelser,
ikke to uafhængige patienter, raske kontroller eller biologiske replikater.
Puromycin/NMD-hæmning kan påvirke hele transkriptomet. Udgangspunktet er
derfor separate, protokolmatchede nt- og puromodeller. En fælles model
med parring og behandling kræver særskilt understøttelse i valgt
softwareversion og dokumenteret sammenlignende validering; dette antages
ikke som en eksisterende funktion.

## Hvor mange prøver og kontroller?

* **Nu:** bevar alle individuelt indeksede BAM'er og metadata. Den
  eksisterende pooled reference kan være nyttig i IGV, men er ikke 20 eller
  et andet antal uafhængige kontroller og kan ikke træne en individbaseret
  outliermodel.
* **Pilot:** planlæg mindst **50 QC-godkendte, uafhængige personer pr.
  homogen modelkohorte**. FRASER-vignetten anbefaler mindst 50 prøver fra
  samme væv og protokol. OUTRIDER-artiklen fandt, at omtrent 50-60 prøver
  var tilstrækkeligt til at genfinde de fleste stærke ekspressionsoutliers
  i netop deres kohorte. Dette er litteraturbaserede pilotmål, ikke et
  universelt minimum eller klinisk acceptgrænse.
* **Større kohorte:** sigt eksempelvis mod 100+ personer pr. betingelse
  og undersøg fortsat forbedring i følsomhed og robusthed. Tallet er et
  planlægningsmål, ikke evidens for tilstrækkelig klinisk performance.
  Ekskludér dårlige prøver og undersøg batch- og behandlingskonfounding;
  et nominelt antal kompenserer ikke for inkompatible prøver.
* **Uafhængig evaluering:** reserver patientuafhængige cases, helst også
  senere batches, der ikke har været brugt til valg af model, tærskel eller
  genfilter. Antal fastlægges ved forudbestemte acceptkriterier, observeret
  performance og konfidensintervaller – ikke ved at kalde de 50
  træningsprøver for validering.
* **Kliniske referencekontroller:** raske protokolmatchede donorer og/eller
  egnede patientprøver uden kendt relevant effekt i det aktuelle gen kan
  være nyttige. Patienter er ikke automatisk raske kontroller. Håndtér
  kendt sygdomsvariant, beslægtning, fælles batches, manglende WGS/WES-
  dækning og mulige uopdagede effekter; fastlæg og dokumentér udvælgelsen.
  Flere referencer øger ikke automatisk evidensstyrken.

Et eksempel på usikkerhed: hvis 59 **uafhængige, repræsentative** kendt
positive cases alle påvises, er den eksakte ensidige 95 % nedre
konfidensgrænse for sensitivitet ca. 95 %. Det er ikke et anbefalet
universelt antal, og det siger intet om specificitet eller hver enkelt
effekttype. Specificitet, ikke-vurderbare analyser, varianter med forskellig
effekt og repræsentativitet skal vurderes separat.

### Aktuelt udgangspunkt og tidsestimat

De cirka 30 er ifølge de oplyste data patienter, og begge betingelser
findes for alle eller næsten alle. De cirka 30 personer har ikke varianter
i de samme gener og kan derfor *potentielt* indgå som genrelevante
referencer for hinanden; de er ikke beskrevet som en separat kohorte af
raske donorer. Det er ikke tilstrækkeligt blot at antage, at et gen er
upåvirket: tjek dækningsområde, variantdata, beslægtning, fænotype,
batch, protokol, QC og eventuel kendt effekt. Antallet af egnede
referencer vil variere pr. gen, og den aktuelle patient må ikke tælles
som sin egen uafhængige kontrol. De 30 patientpar må ikke tælles igen
som 30 ekstra raske kontrolpar.

Hvis alle 30 er uafhængige og består QC i begge betingelser, og der
fremover kommer fire brugbare nye patientpar per måned, er 50 par omtrent
fem måneder væk og 100 par omtrent 17-18 måneder væk. Manglende betingelser,
eksklusioner eller uegnede referencer forlænger tidsplanen. Tidsestimatet
vedrører modelkohorten; uafhængig klinisk validering kommer derudover.

## Fase 1: inventar, metadata og sikkerhed

1. Opret et adgangsbegrænset manifest: pseudonymt patient-ID, prøve-ID,
   nt/puro, væv/celletype, dyrkning, puromycinprotokol og tidspunkt,
   RNA-kvalitet, biblioteksmetode og strandedness, sekventeringsdato,
   flowcelle, læselængde, dybde, batches, aligner og referenceversion.
   Registrér relationer/duplikater og hvilke prøver der er kontroller.
2. Registrér WGS/WES-reference, varianter, kvalitets- og dækningsstatus
   samt transkriptversion. Ingen WES-dækning eller manglende variant er
   ikke det samme som en dokumenteret negativ genomisk kontrol.
3. Gem individuelle BAM/BAI og kvalitetssikrede FASTQ efter lokal
   opbevaringspolitik. Registrér checksums, softwareversioner, parametre,
   GTF/FASTA/cache-checksums og sporbarhed fra analyse til rapport.
   Adgangsstyring, logning og datadeling skal følge lokale regler. Send
   ikke patientvarianter til eksterne tjenester uden godkendelse.
4. Kontroller køn/identitet og match mellem RNA og WGS/WES efter godkendt
   metode. Definér håndtering af prøveforbytning, kontaminering, beslægtning
   og gentagne prøver; undgå at ét individ tæller flere gange i en model.

## Fase 2: harmonisering og QC

1. Fastfrys GRCh38-reference, annotation og én reproducerbar tællemetode
   og softwareversion for hver model. Behandl historiske prøver identisk,
   helst fra FASTQ hvis realignment er nødvendigt; genbrug kun BAM, hvis
   alignerings- og referencekompatibilitet er dokumenteret.
2. Følg pr. prøve mapping, antal informative læsninger, rRNA, duplikater,
   exon/intron-fordeling, gen-/transkriptdækning, strand, junction-counts,
   RNA-kvalitet og 5'/3'-bias. Definér grænser på egne data med kendte
   positive og negative cases; ikke alle metrikker passer til total RNA.
3. Plot PCA, bibliotekstørrelse og kendte batches separat for nt/puro.
   Undersøg, om behandling, dato eller batch er fuldstændig konfounderet
   med patientgruppe; latent korrektion kan ikke garantere at dette løses.
   Foruddefinér eksklusion og håndtering af protokolskift.

## Fase 3: pilot OUTRIDER og FRASER

1. Vælg aktuelle, reproducerbart installerede Bioconductor-versioner, lås
   afhængigheder og opbyg separate modeldatasæt for nt og puro. OUTRIDER
   kræver sammenlignelige gen-counts og filtrering af utilstrækkeligt
   udtrykte gener; FRASER kræver individuelt sammenlignelige
   junction-/intron-counts og filtrering af utilstrækkelig dækning.
   Brug pakkernes dokumenterede workflow for den konkrete version.
2. Start med retrospektiv pilot ved cirka 50-60 kvalificerede personer pr.
   betingelse. Kør prøve- og batchkontroller, kontroller kendte syntetiske
   og biologiske positive cases og sammenlign med manuel IGV. Evaluer
   latent model, filtrering og multiple test-korrektion. Opbevar ikke kun
   signifikante resultater: registrér også gen/hændelser med utilstrækkelig
   dækning eller filtrering.
3. Evaluer recall for kendte exon skips, nye donor-/acceptor-sites,
   pseudoexons, retention og ekspressionsfald separat. Metoderne er ikke
   nødvendigvis følsomme for alle disse effekter. Valider også falske
   positiver, batching, transkriptvalg og konsekvenser af NMD-hæmning.
4. Beskrivelig sammenstilling med WGS/WES: annotér påvirket gen, genomisk
   afstand, transkript og kendt variant. Bevar alle outliers i den
   unbiased analyse; prioritér derefter med genomik og fænotype. Angiv
   WES-blinde regioner og undgå cirkelslutning og dobbelt evidenstælling.

## Fase 4: validering og indførelse

* Lav en forhåndsgodkendt protokol med tilsigtet anvendelse,
  inklusionskriterier, versioner, måleområde og acceptgrænser for
  sensitivitet, specificitet, reproducerbarhed, robusthed, ikke-vurderbare
  gener og svartid. Test kendte positive/negative patientcases på tværs af
  relevante effekter, ekspressionsniveauer, sekventerings- og batchforhold.
  Afklar antal med laboratoriets statistiske og kliniske ansvarlige.
* Hold patienter og så vidt muligt batches adskilt mellem udvikling og
  holdout; alle beslægtede/gentagne prøver for samme individ følger samme
  split. Negative cases skal have passende uafhængig sandhedsreference.
  Kendte positive resultater alene må ikke bruges til at estimere specificitet.
* Dokumentér orthogonal bekræftelse hvor relevant, fx målrettet RT-PCR og
  sekventering, samt begrænsninger ved hvert assay. Opstil rapportskabelon
  med observation, prædiktion, dækningsgrundlag, alternativ forklaring og
  usikkerhed hver for sig. RNA-outlier er ikke automatisk PS3/BS3; undgå
  dobbelt anvendelse af samme evidens.
* Fastfrys og versionsmærk modelkohorte, tællepipeline, annotationer og
  tærskler. En ny patient må ikke umærket påvirke sin egen reference.
  Vurder prospektiv inference ud fra valgt softwareversions faktiske
  understøttelse; hvis kohorten genberegnes sammen med nye prøver,
  dokumentér inklusion, datalækager og effekten af opdateringen.
  Sammenlign gammel/ny model før kontrolleret frigivelse.
* Kør først i skyggeforløb, gennemgå diskordanser og falske fund, og indfør
  først klinisk efter lokal analytisk og klinisk validering, godkendt
  kvalitetssystem, SOP, ansvarlig signering og relevant regulatorisk
  vurdering. Overvåg drift, batch, andel ikke-vurderbare og ændringer over
  tid; definér stopkriterier og tilbagekaldelse/reanalyse.

## Kilder

* OUTRIDER: Brechtmann et al., *Am J Hum Genet* 2018,
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6288422/
* FRASER-vignette, anbefaling om kohorte og metode:
  https://bioconductor.org/packages/release/bioc/vignettes/FRASER/inst/doc/FRASER.pdf
* OUTRIDER-vignette: https://bioconductor.posit.co/packages/3.24/bioc/vignettes/OUTRIDER/inst/doc/OUTRIDER.pdf
* ClinGen SVI Splicing Subgroup, RNA-evidens:
  https://pubmed.ncbi.nlm.nih.gov/37352859/
* Regtools junctionformat og filtre:
  https://regtools.readthedocs.io/en/stable/commands/junctions-extract/
