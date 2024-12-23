# Single-Cell RNA Sequencing (scRNA-seq) Analysis
Overview
Single-cell RNA sequencing (scRNA-seq) is a powerful technique that enables the measurement of gene expression at the level of individual cells, offering a higher resolution view of cellular diversity compared to traditional bulk RNA sequencing. While bulk RNA-seq measures the average gene expression across a population of cells, scRNA-seq isolates and analyzes RNA from each cell separately, allowing researchers to uncover subtle differences in gene expression that are often masked in pooled samples. This ability to capture cellular heterogeneity is particularly valuable for identifying rare cell types, understanding complex developmental processes, and exploring disease mechanisms. The typical workflow involves isolating individual cells, capturing their RNA, converting it into complementary DNA (cDNA), and sequencing it to generate vast amounts of data. Advanced bioinformatics tools are then used to process, analyze, and interpret the resulting data, often revealing previously unrecognized cell types, states, and regulatory networks. Although scRNA-seq provides rich insights into cellular function, it also presents challenges, such as data sparsity, technical noise, and high costs. Nevertheless, its applications are vast, ranging from cancer research and immune system studies to developmental biology and drug discovery, making it a transformative tool in modern genomics.

# Projects and Select Sample Outputs (for real output can be found in notebooks)
## Brief Dataset & Project Description: 
Clustering of Peripheral Blood Mononuclear Cells (PBMCs). - PBMC3k.ipynb
Here we will be analyzing a dataset of Peripheral Blood Mononuclear Cells (PBMC).There are 2,700 single cells that were sequenced on the Illumina NextSeq 500. Peripheral Blood Mononuclear Cells (PBMCs) are a diverse group of immune cells found in peripheral blood, essential for immune surveillance and response. They consist primarily of cells derived from two major lineages: lymphoid and myeloid. The lymphoid components include T cells (involved in adaptive immunity), B cells (responsible for antibody production), and Natural Killer (NK) cells (which target infected or cancerous cells). The myeloid component includes monocytes, which can differentiate into macrophages or dendritic cells to aid in pathogen detection and immune activation. PBMCs are crucial for the body's defense against infections and diseases, and they are widely used in immunological research to study immune function, disease mechanisms, and therapeutic interventions.

## Select Results: 

<img width="346" alt="Screenshot 2024-12-20 at 8 12 48 PM" src="https://github.com/user-attachments/assets/a85818b6-dbb1-4d83-8472-1f33993e935d">

Fig. UMAP of single-cell PBMC data, with cells colored by  and labeled by assigned cell type. Clusters represent distinct cell populations based on gene expression.
## Brief Dataset & Project Description:
## Clustering of Non-small cell lung cancer (NSCLC) dissociated tumor cells - NSCLC.ipynb
Non-small cell lung cancer (NSCLC) dissociated tumor cells from 7 donors were obtained from Discovery Life Sciences. Cells were labeled with TotalSeq™-B Human TBNK Cocktail (BioLegend). Each donor was CellPlexed and pooled at equal proportions. Viable cells in the pool were identified by 7AAD staining and sorted via FACS.
Gene Expression and CellPlex libraries were generated from ~33,000 cells as described in the Chromium Single Cell 3' Reagent Kits User Guide (v3.1 Chemistry Dual Index) with Feature Barcode technology for Cell Surface Protein and Cell Multiplexing (CG000390 Rev B) using the Chromium X and sequenced on an Illumina NovaSeq 6000 to a read depth of approximately 70,000 mean reads per cell for Gene Expression and 25,000 mean reads per cell for CellPlex.

## Select Results: 

<img width="522" alt="Screenshot 2024-12-20 at 10 38 54 PM" src="https://github.com/user-attachments/assets/fc494e74-ff38-465e-94c2-a80455f9249b" />

Fig. UMAP of single-cell Non-small cell lung cancer (NSCLC) dissociated tumor cells data 


