"""
===================================================================================
AMBIENT RNA REMOVAL - DATASET REFERENCE GUIDE
===================================================================================

This file contains code snippets for downloading and loading various scRNA-seq
datasets suitable for testing ambient RNA removal methods.

Datasets are organized by:
1. Built-in datasets (Scanpy/CellRank - easy to download)
2. Public repositories (10X Genomics, GEO, etc.)
3. Recommended for ambient RNA testing

Author: Research meeting prep
Date: October 2025
===================================================================================
"""

import anndata
import numpy as np
import scanpy as sc
from scipy.sparse import issparse

# ===================================================================================
# SECTION 1: BUILT-IN DATASETS (Scanpy/CellRank)
# ===================================================================================

def load_pbmc3k():
    """
    PBMC 3K - Peripheral Blood Mononuclear Cells
    
    Source: 10X Genomics
    Cells: ~3,000
    Platform: 10X Chromium
    
    ⚠️  NOTE: This is PRE-FILTERED. For ambient RNA testing, you need RAW counts!
    Use this for: Quick testing, method demonstration
    Not ideal for: Ambient RNA validation (already filtered)
    """
    print("📥 Loading PBMC3K dataset...")
    adata = sc.datasets.pbmc3k_processed()
    print(f"   Shape: {adata.shape}")
    print(f"   Cell types: {adata.obs['louvain'].nunique()} clusters")
    return adata


def load_pbmc3k_raw():
    """
    PBMC 3K - RAW UNFILTERED VERSION
    
    This loads the raw counts BEFORE filtering - better for ambient RNA testing!
    """
    print("📥 Loading PBMC3K RAW (unfiltered)...")
    adata = sc.datasets.pbmc3k()  # No '_processed' suffix = raw data
    print(f"   Shape: {adata.shape}")
    print(f"   Total UMI range: {adata.X.sum(axis=1).min():.0f} - {adata.X.sum(axis=1).max():.0f}")
    return adata


def load_paul15_bone_marrow():
    """
    Paul et al. 2015 - Mouse Bone Marrow
    
    Source: Paul et al. (2015) Cell
    Cells: ~2,700
    System: Hematopoiesis (blood cell development)
    
    Good for: High contamination scenarios
    Why: Bone marrow has diverse cell types + high ambient RNA
    """
    print("📥 Loading Paul15 bone marrow dataset...")
    adata = sc.datasets.paul15()
    print(f"   Shape: {adata.shape}")
    print(f"   Cell types: {adata.obs['paul15_clusters'].nunique()} clusters")
    return adata


def load_bone_marrow_cellrank():
    """
    Bone Marrow - CellRank Dataset
    
    Source: CellRank package
    Larger bone marrow dataset with velocity information
    
    Installation: pip install cellrank
    """
    try:
        import cellrank as cr
        print("📥 Loading CellRank bone marrow dataset...")
        adata = cr.datasets.bone_marrow(path='datasets/bone_marrow.h5ad')
        print(f"   Shape: {adata.shape}")
        return adata
    except ImportError:
        print("❌ CellRank not installed. Run: pip install cellrank")
        return None


# ===================================================================================
# SECTION 2: REPROGRAMMING / DEVELOPMENTAL DATASETS
# ===================================================================================

def load_schiebinger_reprogramming():
    """
    Schiebinger et al. 2019 - Cell Reprogramming
    
    Source: Schiebinger et al. (2019) Cell
    Cells: ~106,000
    System: Fibroblast to neuron reprogramming
    Time points: Multiple (timecourse)
    
    Download from: https://github.com/AllonKleinLab/Schiebinger-2019-cell
    Or use CellRank: cr.datasets.reprogramming_schiebinger()
    
    Good for: Temporal dynamics + ambient RNA
    """
    print("📥 Loading Schiebinger reprogramming dataset...")
    
    # Option 1: Load from local file
    try:
        adata = anndata.read_h5ad('datasets/reprogramming_schiebinger.h5ad')
        print(f"   ✅ Loaded from local file")
    except:
        # Option 2: Download from CellRank
        try:
            import cellrank as cr
            print("   Downloading from CellRank...")
            adata = cr.datasets.reprogramming_schiebinger(
                path='datasets/reprogramming_schiebinger.h5ad',
                subset_to_serum=False
            )
        except:
            print("   ❌ Could not load. Download manually from:")
            print("   https://github.com/AllonKleinLab/Schiebinger-2019-cell")
            return None
    
    print(f"   Shape: {adata.shape}")
    print(f"   Time points: {adata.obs['day'].unique() if 'day' in adata.obs else 'N/A'}")
    
    return adata


def filter_schiebinger_data(adata, serum=True, day=None):
    """
    Filter Schiebinger dataset by condition and timepoint
    
    Parameters:
    -----------
    adata : AnnData
        Schiebinger dataset
    serum : bool
        Filter for serum condition
    day : str or None
        Specific day to filter (e.g., '0.0', '2.0', etc.)
    
    Returns:
    --------
    Filtered AnnData object
    """
    print("🔍 Filtering Schiebinger data...")
    
    # Filter by serum
    if serum and 'serum' in adata.obs:
        adata_filtered = adata[adata.obs["serum"] == 'True']
        print(f"   After serum filter: {adata_filtered.shape}")
    else:
        adata_filtered = adata
    
    # Filter by day
    if day is not None and 'day' in adata_filtered.obs:
        adata_filtered = adata_filtered[adata_filtered.obs["day"] == str(day)]
        print(f"   After day {day} filter: {adata_filtered.shape}")
    
    return adata_filtered


def load_morris_reprogramming():
    """
    Morris et al. 2019 - Fibroblast Reprogramming
    
    Source: Biddy et al. (2018) Nature
    Cells: ~104,887
    System: Mouse fibroblast → induced endoderm progenitors
    Platform: 10X Chromium + Drop-seq
    Time points: 8 timepoints (days 0-28)
    
    Installation: pip install cellrank
    """
    try:
        import cellrank as cr
        print("📥 Loading Morris reprogramming dataset...")
        adata = cr.datasets.reprogramming_morris(
            path='datasets/reprogramming_morris.h5ad',
            subset='full'
        )
        print(f"   Shape: {adata.shape}")
        print(f"   Time points: {adata.obs['day'].unique() if 'day' in adata.obs else 'N/A'}")
        return adata
    except ImportError:
        print("❌ CellRank not installed. Run: pip install cellrank")
        return None


# ===================================================================================
# SECTION 3: DEVELOPMENTAL BIOLOGY DATASETS
# ===================================================================================

def load_zebrafish_development():
    """
    Farrell et al. 2018 - Zebrafish Embryogenesis
    
    Source: Farrell et al. (2018) Science
    Cells: ~2,434 (axial mesoderm lineage)
    System: Zebrafish embryo development
    Platform: Drop-seq
    Time points: 12 timepoints (3.3-12 hours post-fertilization)
    
    Good for: Developmental trajectories + ambient RNA
    Installation: pip install cellrank
    """
    try:
        import cellrank as cr
        print("📥 Loading zebrafish development dataset...")
        adata = cr.datasets.zebrafish(path='datasets/zebrafish.h5ad')
        print(f"   Shape: {adata.shape}")
        return adata
    except ImportError:
        print("❌ CellRank not installed. Run: pip install cellrank")
        return None


def load_pancreas_development():
    """
    Bastidas-Ponce et al. 2019 - Mouse Pancreas Development
    
    Source: Bastidas-Ponce et al. (2019) Development
    Cells: ~2,531
    System: Murine pancreas at E15.5 (endocrinogenesis)
    Platform: 10X Chromium
    
    Good for: Developmental biology + organ development
    Installation: pip install cellrank
    
    Versions:
    - 'raw': Raw spliced/unspliced counts
    - 'preprocessed': Filtered, normalized, HVG selected
    - 'preprocessed-kernel': With velocity kernel
    """
    try:
        import cellrank as cr
        print("📥 Loading pancreas development dataset...")
        
        # Choose version
        kind = 'preprocessed'  # or 'raw' or 'preprocessed-kernel'
        
        adata = cr.datasets.pancreas(
            path='datasets/endocrinogenesis_day15.5.h5ad',
            kind=kind
        )
        print(f"   Shape: {adata.shape}")
        print(f"   Version: {kind}")
        return adata
    except ImportError:
        print("❌ CellRank not installed. Run: pip install cellrank")
        return None


# ===================================================================================
# SECTION 4: 10X GENOMICS PUBLIC DATASETS (HIGH QUALITY)
# ===================================================================================

def download_10x_dataset_instructions():
    """
    How to download datasets from 10X Genomics website
    
    These are HIGH QUALITY datasets with raw unfiltered matrices!
    Perfect for ambient RNA testing!
    """
    
    instructions = """
    ============================================================================
    DOWNLOADING 10X GENOMICS PUBLIC DATASETS
    ============================================================================
    
    Website: https://www.10xgenomics.com/datasets
    
    RECOMMENDED DATASETS FOR AMBIENT RNA TESTING:
    
    1. 10k PBMCs from a Healthy Donor (v3 chemistry)
       - Cells: ~10,000
       - Link: https://www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0
       - Files needed: "Feature / cell matrix HDF5 (raw)"
       
    2. 10k Brain Cells from an E18 Mouse (v3 chemistry)
       - Cells: ~10,000
       - Link: https://www.10xgenomics.com/datasets/10-k-brain-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0
       - Good for: High contamination (brain has fragile cells)
       
    3. 10k Human PBMCs, 3' v3.1, Chromium X
       - Cells: ~10,000
       - Link: https://www.10xgenomics.com/datasets/10-k-human-pbm-cs-3-v-3-1-chromium-x-3-1-high
       - Latest chemistry
    
    4. 20k Human PBMCs, 3' HT v3.1, Chromium X
       - Cells: ~20,000
       - More cells for robust testing
    
    DOWNLOAD INSTRUCTIONS:
    ----------------------
    1. Click on dataset link
    2. Find "Feature / cell matrix HDF5 (raw)"
    3. Download the .h5 file
    4. Load in Python:
    
        import scanpy as sc
        adata = sc.read_10x_h5('path/to/downloaded_file.h5')
    
    WHY THESE ARE GOOD:
    -------------------
    ✅ RAW unfiltered matrices (includes empty droplets)
    ✅ High quality (10X official datasets)
    ✅ Well-documented
    ✅ Multiple tissue types available
    ✅ Different chemistries for comparison
    
    ============================================================================
    """
    
    print(instructions)
    return instructions


def load_10x_h5_file(filepath):
    """
    Load a 10X Genomics HDF5 file
    
    Parameters:
    -----------
    filepath : str
        Path to .h5 file downloaded from 10X website
    
    Returns:
    --------
    AnnData object with raw counts
    """
    print(f"📥 Loading 10X dataset from {filepath}...")
    try:
        adata = sc.read_10x_h5(filepath)
        print(f"   ✅ Loaded successfully")
        print(f"   Shape: {adata.shape}")
        print(f"   Total UMI per cell: {adata.X.sum(axis=1).mean():.0f} (mean)")
        return adata
    except Exception as e:
        print(f"   ❌ Error loading file: {e}")
        return None


# ===================================================================================
# SECTION 5: RECOMMENDED WORKFLOW
# ===================================================================================

def preprocess_for_ambient_testing(adata, min_genes=200, min_cells=3):
    """
    Minimal preprocessing for ambient RNA testing
    
    ⚠️  IMPORTANT: Don't filter too aggressively!
    We need low-quality cells and empty droplets for ambient estimation!
    
    Parameters:
    -----------
    adata : AnnData
        Raw count matrix
    min_genes : int
        Minimum genes per cell (keep this LOW)
    min_cells : int
        Minimum cells per gene
    """
    print("🔧 Preprocessing for ambient RNA testing...")
    print(f"   Initial shape: {adata.shape}")
    
    # Basic gene filtering only
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print(f"   After gene filter: {adata.shape}")
    
    # Calculate QC metrics but DON'T filter yet
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    print(f"   UMI range: {adata.obs['total_counts'].min():.0f} - {adata.obs['total_counts'].max():.0f}")
    print(f"   Gene range: {adata.obs['n_genes_by_counts'].min():.0f} - {adata.obs['n_genes_by_counts'].max():.0f}")
    
    # Identify cells vs empty droplets (but keep both!)
    adata.obs['is_cell'] = adata.obs['total_counts'] > 500  # Adjust threshold
    print(f"   Cells: {adata.obs['is_cell'].sum()}")
    print(f"   Empty droplets: {(~adata.obs['is_cell']).sum()}")
    
    print("   ✅ Preprocessing complete - kept empty droplets for ambient estimation!")
    
    return adata


# ===================================================================================
# SECTION 6: EXAMPLE USAGE
# ===================================================================================

def example_workflow():
    """
    Example workflow for loading and preparing data for ambient RNA testing
    """
    
    print("\n" + "="*80)
    print("EXAMPLE WORKFLOW: Loading Data for Ambient RNA Testing")
    print("="*80 + "\n")
    
    # Option 1: Built-in dataset (quick test)
    print("OPTION 1: Quick test with built-in data")
    print("-" * 40)
    adata_pbmc = load_pbmc3k_raw()
    if adata_pbmc:
        adata_pbmc = preprocess_for_ambient_testing(adata_pbmc)
    
    print("\n")
    
    # Option 2: Bone marrow (high contamination)
    print("OPTION 2: Bone marrow (high contamination)")
    print("-" * 40)
    adata_bm = load_paul15_bone_marrow()
    
    print("\n")
    
    # Option 3: Download 10X data
    print("OPTION 3: 10X Genomics public datasets")
    print("-" * 40)
    download_10x_dataset_instructions()
    
    print("\n")
    print("="*80)
    print("✅ Examples complete!")
    print("="*80)


# ===================================================================================
# DATASET SUMMARY TABLE
# ===================================================================================

DATASET_SUMMARY = """
╔══════════════════════════════════════════════════════════════════════════════╗
║                      DATASET SUMMARY FOR AMBIENT RNA TESTING                  ║
╚══════════════════════════════════════════════════════════════════════════════╝

┌─────────────────────┬──────────┬────────────────────┬──────────────────────┐
│ Dataset             │ Cells    │ Good For           │ Contamination Level  │
├─────────────────────┼──────────┼────────────────────┼──────────────────────┤
│ PBMC3K (raw)        │ ~3,000   │ Quick testing      │ Low-Medium          │
│ Paul15 Bone Marrow  │ ~2,700   │ Diverse cell types │ Medium-High         │
│ Brain 10k           │ ~10,000  │ Fragile cells      │ High ⭐             │
│ Schiebinger         │ ~106,000 │ Time-series        │ Medium              │
│ Morris              │ ~104,000 │ Reprogramming      │ Medium              │
│ Zebrafish           │ ~2,400   │ Development        │ Medium              │
│ Pancreas            │ ~2,500   │ Organogenesis      │ Medium              │
└─────────────────────┴──────────┴────────────────────┴──────────────────────┘

⭐ = Highly recommended for ambient RNA testing

RECOMMENDED DATASETS (in order of priority):
1. 10X Brain 10k - High contamination, large, well-documented
2. Bone marrow (Paul15 or CellRank) - Diverse cell types
3. Schiebinger - Large sample size, time-series
4. PBMC3K raw - Quick validation

KEY REQUIREMENTS FOR AMBIENT RNA TESTING:
✅ RAW unfiltered counts
✅ Empty droplets included
✅ Diverse cell populations
✅ Known to have contamination issues
"""


def print_dataset_summary():
    """Print the dataset summary table"""
    print(DATASET_SUMMARY)


# ===================================================================================
# MAIN EXECUTION
# ===================================================================================

if __name__ == "__main__":
    print("\n" + "="*80)
    print("🧬 AMBIENT RNA REMOVAL - DATASET REFERENCE GUIDE")
    print("="*80 + "\n")
    
    print_dataset_summary()
    
    print("\n📋 Run example_workflow() to see how to load datasets")
    print("💡 Tip: Use load_10x_h5_file() for custom 10X downloads\n")
    
    # Uncomment to run example
    # example_workflow()
