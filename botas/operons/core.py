# botas/operons/core.py

def infer_operons(
    bam: str,
    gff: str,
    min_score: float = 0.5,
    classifier: str = "rules",
):
    """
    Core operon inference logic.
    This is intentionally CLI-free.
    """
    print(f"[botas] infer_operons called")
    print(f"  BAM: {bam}")
    print(f"  GFF: {gff}")
    print(f"  min_score: {min_score}")
    print(f"  classifier: {classifier}")

    # TODO: real operon logic
    return None
