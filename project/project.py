# DNA/ RNA pairing
# Transcription and reverse transcription
# Translation
# Nucleotide counter / CG content

def main():

    with open("input.txt", "r") as f1:
        seq = f1.read().upper()
        seq = seq.strip()
        seq = str(seq)
        # screen the sequence
        if " " in seq:
            raise ValueError("Gaps in sequence")
        # raise error if the sequence is neither DNA nor RNA
        if all(nuc in "ATCG" for nuc in seq):
            pass
        elif all(nuc in "AUCG" for nuc in seq):
            pass
        else:
            raise ValueError("Invalid sequence")

    options = [
        "1. DNA/ RNA pairing",
        "2. Transcription/ Reverse transcription",
        "3. RNA Translation",
        "4. Nucleotide counter",
        "5. GC content analysis"
        ]
    print("\nChoose the function by typing the function number\n")
    for option in options:
        print(option + "\n")

    functions = {
        "1": pair,
        "2": transcription,
        "3": translation,
        "4": count,
        "5": GC
        }
    mode = input("Function: ")

    if mode in functions:
        newseq = functions[mode](seq)
    else:
        raise ValueError("No such function")

    with open("output.txt", "w") as f2:
        f2.write("Result: " + "\n")
        f2.write(newseq)

    print("\nCompleted, see results in output.txt\n")

# DNA and RNA pairing
def pair(seq):
    # if DNA
    if "T" in seq and "U" not in seq:
        seq = seq.replace("C", "g").replace("G", "c")
        seq = seq.replace("A", "t").replace("T", "a")
        seq = seq.upper()
    # if RNA
    if "U" in seq and "T" not in seq:
        seq = seq.replace("C", "g").replace("G", "c")
        seq = seq.replace("A", "u").replace("U", "a")
        seq = seq.upper()
    return seq

# transcription and reverse transcription
def transcription(seq):
    # if DNA, transcription
    if "T" in seq and "U" not in seq:
        return seq.replace("T", "U")
    # if RNA, reverse transcription
    if "U" in seq and "T" not in seq:
        return seq.replace("U", "T")

# translation
def translation(seq):
    codons = {
        'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V',
        'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V',
        'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
        'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V',
        'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A',
        'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
        'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
        'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
        'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
        'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
        'UAA': '_Stop_', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
        'UAG': '_Stop_', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
        'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G',
        'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
        'UGA': '_Stop_', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
        'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'
        }
    # screen for DNA
    if not all(nuc in "AUCG" for nuc in seq):
        raise ValueError("Invalid RNA sequence")
    protein = ""
    for nuc in range(0, len(seq), 3):
        codon = seq[nuc:nuc+3]
        aa = codons.get(codon,"")
        protein += aa
    return protein

def count(seq):
    #if DNA
    if "T" in seq and "U" not in seq:
        return f"A: {seq.count('A')}, C: {seq.count('C')}, T: {seq.count('T')}, G: {seq.count('G')}"
    #if RNA
    elif "U" in seq and "T" not in seq:
        return f"A: {seq.count('A')}, C: {seq.count('C')}, U: {seq.count('U')}, G: {seq.count('G')}"

def GC(seq):
    GC_cal = 100 * (seq.count("C") + seq.count("G")) / len(seq)
    GC_per = f"GC content: {GC_cal}%"
    return GC_per

if __name__ == "__main__":
    main()
