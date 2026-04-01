# ARES-GBM Project - TODO List

## 📋 Panoramica

Questo documento elenca tutte le attività da completare per il progetto ARES-GBM e la sottomissione del manuscript a Bioinformatics.

**Ultimo aggiornamento:** 2026-04-01

---

## ✅ COMPLETATO

### Repository GitHub
- [x] Inizializzato repository Git
- [x] Configurato remote GitHub
- [x] Push di tutti gli script
- [x] Dati reali TCGA-GBM caricati
- [x] README completo con badge e documentazione

### Script Implementati (10 totali)
- [x] `flux_analysis.py` - Analisi flussi metabolici
- [x] `gene_ko.py` - Simulazione knockout genici
- [x] `survival_analysis.py` - Analisi sopravvivenza Kaplan-Meier
- [x] `drug_sensitivity.py` - Predizione sensibilità farmaci
- [x] `subtype_analysis.py` - Classificazione sottotipi GBM
- [x] `synthetic_lethality.py` - Screen combinazioni sinergiche
- [x] `drug_repositioning.py` - Screen farmaci FDA-approved
- [x] `biomarker_discovery.py` - Identificazione biomarcatori
- [x] `manuscript_figures.py` - Generazione figure manuscript
- [x] `utils.py` - Funzioni utility

### Dati
- [x] `data/tcga_gbm_expression.csv` - 40 TCGA-GBM + 8 GTEx normal
- [x] `data/tcga_gbm_clinical.csv` - Dati clinici con sopravvivenza
- [x] `data/expression_example.csv` - Dati di esempio
- [x] `data/survival_example.csv` - Dati sopravvivena esempio

### Documenti Manuscript
- [x] `manuscript/main.tex` - Manuscript LaTeX completo
- [x] `manuscript/references.bib` - 18 referenze bibliografiche
- [x] `manuscript/cover_letter_bioinformatics.docx.txt` - Cover letter
- [x] `manuscript/author_contributions_crediT.md` - CRediT taxonomy
- [x] `manuscript/data_availability_statement.md` - Data availability
- [x] `manuscript/conflict_of_interest_disclosure.md` - Conflict of interest
- [x] `manuscript/SUPPLEMENTARY_MATERIALS.md` - Supplementary materials
- [x] `manuscript/BIOINFORMATICS_SUBMISSION_CHECKLIST.md` - Checklist

### Docker Support
- [x] `Dockerfile` - Containerizzazione
- [x] `docker-compose.yml` - Multi-container setup
- [x] `.dockerignore` - Docker ignore rules

### Figure Generate
- [x] `results/figure1_benchmark.pdf` - Benchmark comparison
- [x] `results/figure2_validation.pdf` - Validation results
- [x] `results/figure3_pathways.pdf` - Pathway enrichment
- [x] `results/figure4_network.pdf` - Network topology
- [x] `results/figure5_literature.pdf` - Clinical validation

### Risultati Analisi
- [x] `results/drug_sensitivity_demo.csv` - Drug sensitivity predictions
- [x] `results/drug_repositioning.csv` - Drug repositioning candidates
- [x] `results/subtype_assignments.csv` - Subtype classifications
- [x] `results/subtype_differential_analysis.csv` - Differential fluxes
- [x] `results/synthetic_lethality.csv` - Synthetic lethal pairs
- [x] `results/biomarkers_diagnostic.csv` - Diagnostic biomarkers
- [x] `results/biomarkers_prognostic.csv` - Prognostic biomarkers

---

## 🔄 DA FARE - PRIORITÀ ALTA

### 1. bioRxiv Preprint
- [ ] Attendere approvazione bioRxiv (MS ID: BIORXIV/2026/715908)
- [ ] Verificare email di approvazione (24-72 ore)
- [ ] Ottenere DOI da bioRxiv
- [ ] Aggiornare manuscript con bioRxiv DOI

**Stato:** In attesa di approvazione  
**Deadline:** 2026-04-04

---

### 2. Zenodo DOI per il Codice
- [ ] Creare account Zenodo (se non esiste)
- [ ] Collegare repository GitHub a Zenodo
- [ ] Creare release v1.0.0 su GitHub
- [ ] Ottenere DOI da Zenodo
- [ ] Aggiornare README con DOI Zenodo
- [ ] Aggiornare manuscript con DOI Zenodo

**Stato:** Da iniziare  
**Deadline:** 2026-04-08

---

### 3. Completare Manuscript LaTeX
- [ ] Compilare `manuscript/main.tex` in PDF
  - Opzione A: Installare MiKTeX/TeX Live locale
  - Opzione B: Usare Overleaf (online)
- [ ] Verificare che tutte le figure siano referenziate correttamente
- [ ] Contare le parole (target: 4000-5000)
- [ ] Verificare formato referenze (Vancouver style)
- [ ] Controllare che tutti gli autori siano allineati

**Stato:** Da iniziare  
**Deadline:** 2026-04-10

---

### 4. Preparare Supplementary Materials
- [ ] Generare Supplementary Table S1 (differential fluxes completi)
- [ ] Generare Supplementary Table S2 (gene knockout results)
- [ ] Generare Supplementary Table S3 (survival analysis completa)
- [ ] Generare Supplementary Table S4 (literature curation)
- [ ] Preparare Supplementary Figures S1-S6
- [ ] Creare ZIP con tutti i supplementary files

**Stato:** Da iniziare  
**Deadline:** 2026-04-10

---

### 5. Compilare Documenti per Sottomissione
- [ ] Completare campi `[ ]` in `cover_letter_bioinformatics.docx.txt`
  - [ ] Nomi autori
  - [ ] Affiliazioni
  - [ ] Email corrispondente
  - [ ] ORCID IDs
- [ ] Completare `author_contributions_crediT.md`
  - [ ] Firme digitali di tutti gli autori
- [ ] Completare `conflict_of_interest_disclosure.md`
  - [ ] Firme di tutti gli autori
- [ ] Convertire file `.txt` e `.md` in `.docx` e `.pdf`

**Stato:** Da iniziare  
**Deadline:** 2026-04-10

---

## 🟡 PRIORITÀ MEDIA

### 6. Verifica Finale Manuscript
- [ ] Controllo ortografia e grammatica
- [ ] Verificare coerenza nomenclatura (ARES-GBM vs ARES-SENTINEL-GBM)
- [ ] Controllare che tutte le figure siano citate nel testo
- [ ] Verificare che i risultati nel testo corrispondano ai CSV
- [ ] Far revisionare da co-autori

**Stato:** Da iniziare  
**Deadline:** 2026-04-12

---

### 7. Preparare Figure ad Alta Risoluzione
- [ ] Esportare tutte le figure in TIFF 300+ DPI
- [ ] Verificare dimensioni (larghezza: 8.3cm, 12cm, 17.8cm)
- [ ] Controllare che i font siano leggibili
- [ ] Salvare in cartella `manuscript/figures/`

**Stato:** Da iniziare  
**Deadline:** 2026-04-12

---

### 8. Suggested Reviewers
- [ ] Identificare 3-5 revisori potenziali
- [ ] Verificare conflitti di interesse
- [ ] Preparare lista con email e affiliazioni
- [ ] Inserire nel sistema di submission

**Stato:** Da iniziare  
**Deadline:** 2026-04-15

---

## 🟢 PRIORITÀ BASSA

### 9. Ottimizzazioni Script (Opzionale)
- [ ] Aggiungere test unitari (`tests/`)
- [ ] Implementare logging invece di print
- [ ] Aggiungere type hints completi
- [ ] Creare documentazione Sphinx
- [ ] Configurare CI/CD con GitHub Actions

**Stato:** Opzionale  
**Deadline:** Dopo sottomissione

---

### 10. Materiali Aggiuntivi (Opzionale)
- [ ] Creare video tutorial per YouTube
- [ ] Preparare demo interattiva su Google Colab
- [ ] Creare sito web per il progetto
- [ ] Preparare thread Twitter per promozione

**Stato:** Opzionale  
**Deadline:** Dopo accettazione

---

## 📅 Sottomissione a Bioinformatics

### Checklist Pre-Submission
- [ ] Manuscript PDF pronto
- [ ] Supplementary materials pronti
- [ ] Cover letter completata e firmata
- [ ] Author contributions compilate
- [ ] Conflict of interest compilato
- [ ] Data availability statement pronto
- [ ] Figure ad alta risoluzione pronte
- [ ] bioRxiv DOI ottenuto
- [ ] Zenodo DOI ottenuto
- [ ] Reviewers suggeriti identificati

### Processo di Submission
1. [ ] Andare su https://mc.manuscriptcentral.com/bioinformatics
2. [ ] Creare account (se non esiste)
3. [ ] Iniziare nuova submission
4. [ ] Selezionare "Original Paper"
5. [ ] Caricare manuscript
6. [ ] Caricare supplementary materials
7. [ ] Caricare figure
8. [ ] Inserire metadata (titolo, abstract, keywords)
9. [ ] Inserire autori e affiliazioni
10. [ ] Inserire suggested reviewers
11. [ ] Caricare cover letter
12. [ ] Revisionare e submit

**Target Submission Date:** 2026-04-15

---

## 📊 Timeline

| Data | Milestone |
|------|-----------|
| 2026-04-04 | Approvazione bioRxiv |
| 2026-04-08 | Zenodo DOI ottenuto |
| 2026-04-10 | Manuscript completo |
| 2026-04-12 | Figure e supplementary pronti |
| 2026-04-15 | **Sottomissione a Bioinformatics** |
| 2026-04-22 | Initial screening completato |
| 2026-05-20 | Peer review responses (stimato) |
| 2026-06-15 | Revisione submit (se necessaria) |
| 2026-07-01 | Accettazione (stimata) |

---

## 🔗 Link Utili

### Repository e Dati
- **GitHub:** https://github.com/ARES-Sentinel-GBM/-gbm-research
- **bioRxiv:** https://www.biorxiv.org/
- **Zenodo:** https://zenodo.org/
- **TCGA:** https://portal.gdc.cancer.gov/

### Journal
- **Bioinformatics:** https://academic.oup.com/bioinformatics
- **Submission:** https://mc.manuscriptcentral.com/bioinformatics
- **Author Guidelines:** https://academic.oup.com/bioinformatics/pages/General_Instructions

### Strumenti
- **Overleaf:** https://www.overleaf.com/ (LaTeX online)
- **MiKTeX:** https://miktex.org/ (LaTeX Windows)
- **TeX Live:** https://tug.org/texlive/ (LaTeX cross-platform)

---

## 📧 Contatti Importanti

### Journal Office
- **Email:** bioinformatics@oxfordjournals.org

### Supporto Tecnico
- **GitHub Issues:** https://github.com/ARES-Sentinel-GBM/-gbm-research/issues
- **Email progetto:** admin@ares-bio.com

---

## 📝 Note

- Tutti i deadline sono indicativi e possono essere aggiustati
- Priorità può cambiare in base a feedback dei co-autori
- Documentare ogni cambiamento nel repository
- Mantenere backup locale di tutti i file

---

**Ultima revisione:** 2026-04-01  
**Prossima revisione:** 2026-04-08
