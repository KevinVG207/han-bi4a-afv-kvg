from flask import Flask, render_template, request

app = Flask(__name__)


@app.route('/')
def index():
    return render_template("index.html")


@app.route('/dna2prot', methods=["POST", "GET"])
def dna2prot():
    nucleotide_seq = request.args.get("seq", "").lower().replace("u", "t")
    if nucleotide_seq == "":
        protein_seq = ""
    else:
        if not is_dna(nucleotide_seq):
            protein_seq = "ERROR: Not a valid DNA/RNA sequence"
        else:
        # IS valid
            protein_seq = translate(nucleotide_seq)
    return render_template("dna2prot.html", prot=protein_seq)


def is_dna(seq):
    """Deze functie bepaalt of de sequentie (een element uit seqs)
    DNA is.
    Indien ja, return True
    Zo niet, return False
    """
    a_count = seq.count("A")  # Tel alle letters en sla apart op in variabelen
    g_count = seq.count("G")
    c_count = seq.count("C")
    t_count = seq.count("T")
    n_count = seq.count("N")

    total = a_count + g_count + c_count + t_count + n_count  # Bereken totaal lengte van potentieel mRNA

    if total == len(seq):  # Wanneer het totaal overeen komt, is het DNA/mRNA
        return True
    else:
        return False


def translate(seq):
    """
    Neemt een DNA sequentie en transleert deze naar een eiwitsequentie.
    :param seq:
    :return: String met de eiwitsequentie.
    """
    code = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
            'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
            'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
            'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
            'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
            'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
            'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
            'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
            'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
            'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
            'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
            'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
            'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
            'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
            'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
            'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
            }

    seq = seq.lower()
    protein = ""
    try:
        start = seq.index("atg")
        if start != -1:
            seq = seq[start:]
            while len(seq) > 0:
                if len(seq[:3]) < 3:  # Remainder sequence is less than 3 in length
                    break
                curamino = code[seq[:3]]
                if curamino != "*":
                    protein += curamino
                    seq = seq[3:]
                else:
                    protein += curamino
                    break
        return protein
    except ValueError:
        return "ERROR: No start codon found."


if __name__ == '__main__':
    app.run()
