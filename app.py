from flask import Flask, render_template, request
import mysql.connector
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

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


@app.route('/ensembl', methods=["POST", "GET"])
def ensembl():
    if request.method == "POST" and request.form.get("search", "") != "":
        # POST request
        rows = []
        search = request.form.get("search")
        try:
            connection = mysql.connector.connect(host='ensembldb.ensembl.org',
                                                 port=3306,
                                                 user='anonymous',
                                                 db='homo_sapiens_core_95_38')
            cursor = connection.cursor()
            query = """select * from gene where description like '%""" + search + """%' limit 100;"""
            cursor.execute(query)
            rows = cursor.fetchall()
            cursor.close()
            connection.close()
        except mysql.connector.Error as error:
            rows = [error]
        return render_template("ensembl.html", rows=rows)
    else:
        return render_template("ensembl.html", rows=[])


@app.route('/biopython', methods=["POST", "GET"])
def biopython():
    if request.method == "POST" and request.form.get("entry", "") != "":
        # POST request
        url = ""
        user_input = request.form.get("entry").upper().replace("\n", "").replace(" ", "")
        result = ["Input sequence: " + str(user_input)]
        if is_dna(user_input):
            # DNA
            seq_obj = Seq(user_input, IUPAC.ambiguous_dna)
            result += ["DNA Sequence.",
                      "RNA: " + str(seq_obj.transcribe()),
                      "Protein: " + str(seq_obj.translate())]
        elif is_rna(user_input):
            # RNA
            seq_obj = Seq(user_input, IUPAC.ambiguous_rna)
            result += ["RNA Sequence.",
                      "Protein: " + str(seq_obj.translate())]
        elif is_protein(user_input):
            # Protein
            result += ["Protein sequence.",
                      "BLAST to find out the most likely source gene:"]
            url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&QUERY=" + str(user_input)
            # What happens here is that I'm linking to the NCBI BlastP webpage and pre-filling the query box
            # with the given sequence. The user only has to click the "BLAST" button.
            # Doing blasts through the website is much, much faster and means that the webpage
            # won't hang for 15 minutes while the blast happens.
        else:
            # Nothing
            result += ["Not DNA/RNA/Protein"]
        return render_template("biopython.html", result=result, url=url)
    else:
        return render_template("biopython.html", result=[""])


def is_dna(seq):
    """Deze functie bepaalt of de sequentie DNA is.
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


def is_rna(seq):
    """Deze functie bepaalt of de sequentie RNA is.
    Indien ja, return True
    Zo niet, return False
    """
    a_count = seq.count("A")
    g_count = seq.count("G")
    c_count = seq.count("C")
    u_count = seq.count("U")
    n_count = seq.count("N")

    total = a_count + g_count + c_count + u_count + n_count

    if total == len(seq):
        return True
    else:
        return False


def is_protein(seq):
    """Deze functie bepaalt of de sequentie eiwit is.
    Indien ja, return True
    Zo niet, return False
    """
    amino_list = ["A",
                  "R",
                  "N",
                  "D",
                  "B",
                  "C",
                  "E",
                  "Q",
                  "Z",
                  "G",
                  "H",
                  "I",
                  "L",
                  "K",
                  "M",
                  "F",
                  "P",
                  "S",
                  "T",
                  "W",
                  "Y",
                  "V",
                  "X"]
    total = 0
    for amino in amino_list:
        total += seq.count(amino)
    if total == len(seq):
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
