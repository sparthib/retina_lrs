---
config:
  theme: base
  themeVariables:
    fontSize: 22px
  layout: dagre
---
flowchart TD
    %% FASTQ Processing Section
    subgraph FASTQ_Processing["FASTQ Processing"]
        A1[(Raw FASTQ)] --> A2{MinIONQC<br/>Analysis 1}
        A2 --> A3[/YAML Summary & Plots/]
        A1 --> A4{Nanofilt QC<br/>Analysis 2}
        A4 --> A5[(Processed FASTQ)]
        A3 --> A6{MinIONQC Boxplots<br/>Analysis 3}
        A6 --> A7[/Boxplots/]
    end

    %% FASTQ to BAM Section
    subgraph FASTQ_to_BAM["FASTQ to BAM"]
        A5 --> B1{Genome Alignment<br/>Analysis 4}
        A5 --> B2{Transcriptome Alignment<br/>Analysis 5}
        B1 --> B3[(BAM files)]
        B2 --> B4[(BAM files)]
        B3 --> B5{High Quality BAM<br/>Analysis 6}
        B5 --> B6[(Filtered BAM<br/>MAPQ>30)]
        B5 --> B7[/Flagstat Summary/]
    end

    %% BAM QC Section
    subgraph BAM_QC["BAM QC Visualization"]
        B6 --> C1{Flagstat Analysis<br/>Analysis 7}
        C1 --> C1a[/Flagstat Stats/]
        B4 --> C2{Exon-Exon Junction<br/>Analysis 8}
        C2 --> C2a[/Junction Plots/]
        C1a --> C3{Read Type Percentages<br/>Analysis 9 & 10}
        C3 --> C3a[/Read Type Plots/]
    end

    %% Quantification Section
    subgraph Quantification["Quantification"]
        B6 --> D1{Sample-wise Read Class<br/>Analysis 11}
        D1 --> D2[(RDS files)]
        D2 --> D3{Bambu Generate RCS<br/>Analysis 12}
        D3 --> D4[(Extended Annotation GTF<br/>Counts Matrix)]
        D4 --> D5{SQANTI3<br/>Analysis 13}
        D5 --> D6[(Classification.txt<br/>CDS GTF)]
        D2 --> D7{Bambu Quantification<br/>Analysis 14}
        D7 --> D7a[(SE Object)]
        
        B6 --> E1{IsoQuant<br/>Analysis 16}
        E1 --> E2[(IsoQuant Output<br/>Counts & Annotation)]
        
        D4 --> F1{Compare GTFs<br/>Analysis 17}
        E2 --> F1
        F1 --> F2[/Common Isoforms txt/]
        F2 --> D8{Gene Names for Novel<br/>Analysis 15}
        D8 --> D8a[/Gene Names tsv/]
    end

    %% Cleaning Counts Matrix Section
    subgraph Cleaning["Cleaning Counts Matrix"]
        D4 --> G1{Clean Column/Row Names<br/>Analysis 18}
        G1 --> G2[(FT vs RGC & ROs Matrices)]
        G2 --> G3{Filter Common Isoforms<br/>Analysis 19}
        F2 --> G3
        G3 --> G3a[(Filtered Matrices)]
        G3a --> G4{Filter by Gene Biotypes<br/>Analysis 20}
        G4 --> G5[(PCG Gene/Isoform<br/>Counts & CPM)]
    end

    %% DTU DGE DTE Analysis Section
    subgraph DTU_DGE_DTE["DTU DGE DTE Analysis"]
        G5 --> H1{DTE & DGE Analysis<br/>Analysis 21}
        H1 --> H2[/DGE/DTE Results tsv/]
        G5 --> H3{DTU Analysis<br/>Analysis 22}
        D6 --> H3
        H3 --> H4[/DTU Results tsv<br/>Switch Analysis/]
        H4 --> H5{External Protein Analysis<br/>Analysis 23}
        H5 --> H6[/Pfam/SignalP/CPC2 tsv/]
        H2 --> H7{Create Summary<br/>Analysis 24}
        H4 --> H7
        H7 --> H8[/Merged DGE/DTE/DTU tsv/]
    end

    %% Visualization Section
    subgraph Visualization["Visualization"]
        G5 --> I1{PCA Analysis<br/>Analysis 25}
        I1 --> I1a[/PCA Plots/]
        G5 --> I2{Heatmaps<br/>Analysis 26}
        H8 --> I2
        I2 --> I2a[/Heatmap PDFs/]
        H8 --> I3{Volcano Plots<br/>Analysis 27}
        I3 --> I3a[/Volcano Plots/]
        H8 --> I4{RetNet DTU Genes<br/>Analysis 28}
        I4 --> I4a[/RetNet Heatmaps/]
        H4 --> I5{RetNet Switchplots<br/>Analysis 29}
        I4 --> I5
        I5 --> I5a[/RetNet Switchplots/]
        H8 --> I6{Splicing Factor Analysis<br/>Analysis 30}
        I6 --> I6a[/SF Heatmaps/]
        H8 --> I7{Splicing Factor Volcano<br/>Analysis 31}
        I7 --> I7a[/SF Volcano Plots/]
        H4 --> I8{Splicing Factor Switchplots<br/>Analysis 32}
        I6 --> I8
        I8 --> I8a[/SF Switchplots/]
        H8 --> I9{GO Analysis<br/>Analysis 33}
        I9 --> I9a[/GO Dotplots/]
        H8 --> I10{Upset Plots<br/>Analysis 34}
        I10 --> I10a[/Upset Plots/]
        H8 --> I11{DTU-only Upset<br/>Analysis 35}
        I11 --> I11a[/DTU Upset Plot/]
        I11 --> I12{DTU-only GO<br/>Analysis 36}
        I12 --> I12a[/DTU GO Dotplot/]
        G5 --> I13{Short Read Comparison<br/>Analysis 37}
        I13 --> I13a[/Correlation Heatmap/]
        G5 --> I14{Isoforms per Gene<br/>Analysis 38}
        D4 --> I14
        I14 --> I14a[/Isoform Barplot/]
    end

    %% Style definitions
    classDef processing fill:#2196F3,stroke:#0D47A1,stroke-width:2px,color:#fff
    classDef quantification fill:#9C27B0,stroke:#4A148C,stroke-width:2px,color:#fff
    classDef analysis fill:#FF9800,stroke:#E65100,stroke-width:2px,color:#fff
    classDef visualization fill:#4CAF50,stroke:#1B5E20,stroke-width:2px,color:#fff
    classDef dataFile fill:#F44336,stroke:#B71C1C,stroke-width:2px,color:#fff
    classDef plotFile fill:#FFC107,stroke:#F57C00,stroke-width:2px,color:#000

    class FASTQ_Processing,BAM_QC processing
    class Quantification,Cleaning quantification
    class DTU_DGE_DTE analysis
    class Visualization visualization
    class A1,A5,B3,B4,B6,D2,D4,D6,D7a,E2,G2,G3a,G5 dataFile
    class A3,A7,B7,C1a,C2a,C3a,D8a,F2,H2,H4,H6,H8,I1a,I2a,I3a,I4a,I5a,I6a,I7a,I8a,I9a,I10a,I11a,I12a,I13a,I14a plotFile