---
config:
  theme: default
  themeVariables:
    fontSize: 22px
  layout: dagre
---

flowchart TD
    %% Legend
    subgraph Legend["Legend - Data Types (Shapes)"]
        L1[Analysis Step]:::dataAnalysis
        L2([Input Data]):::dataInput
        L3[(Output Data)]:::dataOutput
        L4{{Plot/Visualization}}:::dataPlot
        L5[[Intermediate File]]:::dataIntermediate
        L6{Software Tool}:::dataTool
    end
    
    %% FASTQ Processing Section
    subgraph FASTQ_Processing["FASTQ Processing"]
        style FASTQ_Processing fill:#e1f5fe,stroke:#01579b,stroke-width:3px
        
        A1([Raw FASTQ]):::dataInput
        A2[MinIONQC<br/>Analysis 1]:::dataAnalysis
        A3[[YAML Summary]]:::dataIntermediate
        A3b{{Read Length Plots<br/>Quality Plots}}:::dataPlot
        
        A1 --> A2
        A2 --> A3
        A2 --> A3b
        
        A1 --> A4[Nanofilt QC<br/>Analysis 2]:::dataAnalysis
        A4 --> A5[[Processed FASTQ]]:::dataIntermediate
        
        A3 --> A6[MinIONQC Boxplots<br/>Analysis 3]:::dataAnalysis
        A6 --> A7{{QC Boxplots}}:::dataPlot
    end
    
    %% FASTQ to BAM Section
    subgraph FASTQ_to_BAM["FASTQ to BAM - Alignment"]
        style FASTQ_to_BAM fill:#f3e5f5,stroke:#4a148c,stroke-width:3px
        
        A5 --> B1{Minimap2 Genome<br/>Analysis 4}:::dataTool
        B1 --> B2[[Genome BAM]]:::dataIntermediate
        
        A5 --> B3{Minimap2 Transcriptome<br/>Analysis 5}:::dataTool
        B3 --> B4[[Transcriptome BAM]]:::dataIntermediate
        
        B2 --> B5[High Quality Filter<br/>Analysis 6]:::dataAnalysis
        B4 --> B5
        B5 --> B6[(Filtered BAM<br/>MAPQ≥30)]:::dataOutput
        B5 --> B7[[Flagstat Summary]]:::dataIntermediate
    end
    
    %% BAM QC Section
    subgraph BAM_QC["BAM QC Visualization"]
        style BAM_QC fill:#fff3e0,stroke:#e65100,stroke-width:3px
        
        B6 --> C1[Exon-Exon Junction<br/>Analysis 7]:::dataAnalysis
        C1 --> C2{{Junction Plots}}:::dataPlot
        
        B7 --> C3[Read Type Analysis<br/>Analysis 8-9]:::dataAnalysis
        C3 --> C4{{Read Type Plots}}:::dataPlot
    end
    
    %% Quantification Section
    subgraph Quantification["Quantification"]
        style Quantification fill:#e8f5e9,stroke:#1b5e20,stroke-width:3px
        
        B6 --> D1[Sample-wise Read Class<br/>Analysis 10]:::dataAnalysis
        D1 --> D2[[RDS Files]]:::dataIntermediate
        
        D2 --> D3[Bambu Generate RCS<br/>Analysis 11]:::dataAnalysis
        D3 --> D4[(SE Object<br/>Counts Matrix<br/>Extended GTF)]:::dataOutput
        
        D4 --> D5{SQANTI3<br/>Analysis 12}:::dataTool
        D5 --> D6[(Classification<br/>CDS GTF<br/>Corrected GTF/FASTA)]:::dataOutput
        
        B6 --> D9{IsoQuant<br/>Analysis 14}:::dataTool
        D9 --> D10[(IsoQuant Output<br/>Counts & GTF)]:::dataOutput
    end
    
    %% GTF Comparison Section
    subgraph GTF_Comparison["GTF Comparison"]
        style GTF_Comparison fill:#fce4ec,stroke:#880e4f,stroke-width:3px
        
        D4 --> E1{gffcompare<br/>Analysis 15}:::dataTool
        D10 --> E1
        E1 --> E2[(Common Isoforms)]:::dataOutput
        
        E2 --> E3[Gene Names<br/>Analysis 13]:::dataAnalysis
        E3 --> E4[(Gene Names TSV)]:::dataOutput
    end
    
    %% Counts Matrix Cleaning Section
    subgraph Counts_Cleaning["Counts Matrix Cleaning"]
        style Counts_Cleaning fill:#fce4ec,stroke:#880e4f,stroke-width:3px
        
        D4 --> F1[Clean Counts Matrix<br/>Analysis 16]:::dataAnalysis
        F1 --> F2[(FT vs RGC Matrix<br/>ROs Matrix)]:::dataOutput
        
        F2 --> F3[Filter Common Isoforms<br/>Analysis 17]:::dataAnalysis
        E2 --> F3
        F3 --> F4[[Filtered Matrix]]:::dataIntermediate
        
        F4 --> F5[Filter by Gene Biotype<br/>Analysis 18]:::dataAnalysis
        F5 --> F6[(PCG Counts & CPM)]:::dataOutput
    end
    
    %% DTU/DGE/DTE Analysis Section
    subgraph DTU_DGE_DTE["DTU/DGE/DTE Analysis"]
        style DTU_DGE_DTE fill:#fce4ec,stroke:#880e4f,stroke-width:3px
        
        F6 --> G1[DTE/DGE Analysis<br/>Analysis 19]:::dataAnalysis
        G1 --> G2[(DTE/DGE Results)]:::dataOutput
        
        F6 --> G3[DTU Analysis<br/>Analysis 20]:::dataAnalysis
        D6 --> G3
        G3 --> G4[(DTU Results)]:::dataOutput
        G3 --> G4b{{Switch Plots}}:::dataPlot
        
        G4 --> G5{Protein Analysis<br/>Analysis 21}:::dataTool
        G5 --> G6[(Pfam/SignalP/CPC2)]:::dataOutput
        G6 --> G3
        
        G2 --> G7[Merge Results<br/>Analysis 22]:::dataAnalysis
        G4 --> G7
        G7 --> G8[(Merged DGE/DTE/DTU)]:::dataOutput
    end
    
    %% Visualization Section
    subgraph Visualization["Visualization"]
        style Visualization fill:#fff8e1,stroke:#f57f17,stroke-width:3px
        
        F6 --> H1[PCA Analysis<br/>Analysis 23]:::dataAnalysis
        H1 --> H2{{PCA Plots}}:::dataPlot
        
        F6 --> H3[Heatmaps<br/>Analysis 24]:::dataAnalysis
        G8 --> H3
        H3 --> H4{{Expression Heatmaps}}:::dataPlot
        
        G8 --> H5[Volcano Plots<br/>Analysis 24]:::dataAnalysis
        H5 --> H6{{Volcano Plots}}:::dataPlot
        
        G8 --> H7[RetNet Analysis<br/>Analysis 25]:::dataAnalysis
        H7 --> H8{{IRD Gene Heatmaps}}:::dataPlot
        
        G4 --> H9[RetNet Switch Plots<br/>Analysis 26]:::dataAnalysis
        H7 --> H9
        H9 --> H10{{IRD Switch Plots}}:::dataPlot
        
        G8 --> H11[Splicing Factor Analysis<br/>Analysis 27-29]:::dataAnalysis
        H11 --> H12{{SF Heatmaps & Plots}}:::dataPlot
        
        G8 --> H13[GO Analysis<br/>Analysis 30]:::dataAnalysis
        H13 --> H14{{GO Dotplots}}:::dataPlot
        
        G8 --> H15[Upset Plots<br/>Analysis 31-33]:::dataAnalysis
        H15 --> H16{{Upset & GO Plots}}:::dataPlot
        
        F6 --> H17[Short Read Comparison<br/>Analysis 34]:::dataAnalysis
        H17 --> H18{{Correlation Heatmap}}:::dataPlot
        
        D4 --> H19[Isoforms per Gene<br/>Analysis 35]:::dataAnalysis
        F6 --> H19
        H19 --> H20{{Isoform Distribution}}:::dataPlot
    end
    
    %% Allele Specific Analysis Section
    subgraph ASE["Allele Specific Analysis"]
        style ASE fill:#e3f2fd,stroke:#0d47a1,stroke-width:3px
        
        I1([H9 WGS FASTQ]):::dataInput
        I1 --> I2{Bowtie2<br/>Analysis 36}:::dataTool
        I2 --> I3[[WGS BAM]]:::dataIntermediate
        
        I3 --> I4[Filter BAMs<br/>Analysis 37]:::dataAnalysis
        I4 --> I5[[Filtered WGS BAM]]:::dataIntermediate
        
        I5 --> I6{GATK<br/>Analysis 38}:::dataTool
        I6 --> I7[(Joint VCF)]:::dataOutput
        
        I7 --> I8[VCF Stats<br/>Analysis 39]:::dataAnalysis
        I8 --> I9{{VCF Stats Plots}}:::dataPlot
        
        I7 --> I10{WhatsHap Phase<br/>Analysis 40}:::dataTool
        B6 --> I10
        I10 --> I11[(Phased VCF)]:::dataOutput
        
        I11 --> I12{WhatsHap Haplotag<br/>Analysis 41}:::dataTool
        B6 --> I12
        I12 --> I13[(Haplotagged BAMs<br/>HP1 & HP2)]:::dataOutput
        
        I13 --> I14{featureCounts<br/>Analysis 42}:::dataTool
        I11 --> I14
        I14 --> I15[(ASE Counts Matrix)]:::dataOutput
        
        I15 --> I16[ASE Modeling<br/>Analysis 43]:::dataAnalysis
        I16 --> I17[(ASE DEGs)]:::dataOutput
        I16 --> I17b{{Volcano & GO Plots}}:::dataPlot
        
        I17 --> I18[Chr Distribution<br/>Analysis 44]:::dataAnalysis
        I18 --> I19{{Chr Barplot}}:::dataPlot
        
        B6 --> I20[Variants per Read<br/>Analysis 45]:::dataAnalysis
        I11 --> I20
        I20 --> I21[(Variant Count TSV)]:::dataOutput
        
        I21 --> I22[Variant Percentage<br/>Analysis 46]:::dataAnalysis
        I22 --> I23{{Variant Histograms}}:::dataPlot
    end
    
    %% Connections between main pipeline and ASE
    B6 -.-> I10
    B6 -.-> I12
    B6 -.-> I20
    
    %% Data type styling (consistent colors for all data types)
    classDef dataAnalysis fill:#4CAF50,stroke:#2E7D32,stroke-width:2px,color:#fff
    classDef dataInput fill:#2196F3,stroke:#1565C0,stroke-width:3px,color:#fff
    classDef dataOutput fill:#FF9800,stroke:#E65100,stroke-width:3px,color:#fff
    classDef dataPlot fill:#E91E63,stroke:#AD1457,stroke-width:2px,color:#fff
    classDef dataIntermediate fill:#FFC107,stroke:#F57C00,stroke-width:2px
    classDef dataTool fill:#9C27B0,stroke:#6A1B9A,stroke-width:2px,color:#fff
    
    %% Legend styling
    style Legend fill:#f5f5f5,stroke:#616161,stroke-width:2px