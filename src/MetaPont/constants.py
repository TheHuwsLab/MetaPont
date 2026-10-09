"""Shared schemas and constants for MetaPont."""
from . import MetaPont_Version

UNKNOWN = "unknown"
LINEAGE_RANKS = ("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
LINEAGE_PREFIXES = {
    "Domain": "d__",
    "Kingdom": "k__",
    "Phylum": "p__",
    "Class": "c__",
    "Order": "o__",
    "Family": "f__",
    "Genus": "g__",
    "Species": "s__",
}

# eggNOG-mapper 2.x annotation columns mapped to MetaPont category labels.
FEATURE_COLUMNS = {
    "COG": "COG_category",
    "eggNOG_OG": "eggNOG_OGs",
    "GO": "GOs",
    "EC": "EC",
    "KEGG_KO": "KEGG_ko",
    "KEGG_Pathway": "KEGG_Pathway",
    "KEGG_Module": "KEGG_Module",
    "KEGG_Reaction": "KEGG_Reaction",
    "KEGG_rclass": "KEGG_rclass",
    "BRITE": "BRITE",
    "KEGG_TC": "KEGG_TC",
    "CAZy": "CAZy",
    "BiGG_Reaction": "BiGG_Reaction",
    "PFAM": "PFAMs",
}

EMAPPER_FALLBACK_COLUMNS = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
    "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
    "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
    "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
]

CONTIG_TABLE_COLUMNS = [
    "Sample", "Contig", "Contig_Length", "Mapped_Reads", "Total_Reads",
    "Mapping_Proportion", "Kraken_Classification", "TaxID", "Scientific_Name",
    *LINEAGE_RANKS, "Lineage", "Gene_Count", "Annotated_Gene_Count", "Gene_Read_Support",
    *FEATURE_COLUMNS.keys(),
]
