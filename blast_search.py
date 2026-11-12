import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

# 1. Recorded Search Parameters
QUERY_FILE = "query.fasta"
PROGRAM = "blastn"
DATABASE = "nt"
EXPECT_THRESHOLD = 1e-5
MAX_HITS = 5


def read_fasta(filepath):
    header = ""
    seq = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
            elif line:
                seq.append(line)
    return header, "".join(seq)


def submit_blast(sequence):
    url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    params = {
        "CMD": "Put",
        "PROGRAM": PROGRAM,
        "DATABASE": DATABASE,
        "QUERY": sequence,
        "EXPECT": str(EXPECT_THRESHOLD),
        "HITLIST_SIZE": str(MAX_HITS),
        "FORMAT_TYPE": "XML",
    }
    data = urllib.parse.urlencode(params).encode("utf-8")
    req = urllib.request.Request(url, data=data)
    with urllib.request.urlopen(req) as resp:
        text = resp.read().decode("utf-8")

    rid = None
    rtoe = 15
    for line in text.splitlines():
        if "RID =" in line:
            rid = line.split("=")[1].strip()
        elif "RTOE =" in line:
            rtoe = int(line.split("=")[1].strip())
    return rid, rtoe


def fetch_results(rid, wait_time):
    print(f"Submitted RID: {rid}. Waiting {wait_time}s for NCBI BLAST...")
    time.sleep(wait_time)
    url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID={rid}"

    while True:
        with urllib.request.urlopen(url) as resp:
            xml_text = resp.read().decode("utf-8")
        if "Status=WAITING" in xml_text:
            print("Still running... checking again in 10s.")
            time.sleep(10)
        else:
            return xml_text


def parse_and_save(xml_text, query_header, query_len):
    with open("blast_output.xml", "w", encoding="utf-8") as f:
        f.write(xml_text)

    root = ET.fromstring(xml_text)
    hits = root.findall(".//Hit")

    lines = [
        "# BLAST Search Output & Annotation Table\n",
        "## Search Parameters",
        f"- **Query**: `{query_header}` ({query_len} bp)",
        f"- **Program**: `{PROGRAM}`",
        f"- **Database**: `{DATABASE}` (NCBI Nucleotide)",
        f"- **E-value Threshold**: `{EXPECT_THRESHOLD}`",
        f"- **Max Target Sequences**: `{MAX_HITS}`\n",
        "## Annotation Table of Best Matches\n",
        "| Rank | Accession | Organism / Description | Query Cover | Identity (%) | E-value | Bit Score |",
        "| :--- | :--- | :--- | :---: | :---: | :---: | :---: |",
    ]

    for hit in hits:
        num = hit.findtext("Hit_num")
        acc = hit.findtext("Hit_accession")
        desc = hit.findtext("Hit_def").split(";")[0][:60]
        hsp = hit.find(".//Hsp")

        identity = int(hsp.findtext("Hsp_identity"))
        align_len = int(hsp.findtext("Hsp_align-len"))
        q_from = int(hsp.findtext("Hsp_query-from"))
        q_to = int(hsp.findtext("Hsp_query-to"))
        evalue = float(hsp.findtext("Hsp_evalue"))
        bit_score = float(hsp.findtext("Hsp_bit-score"))

        pct_id = (identity / align_len) * 100
        q_cov = ((abs(q_to - q_from) + 1) / query_len) * 100

        lines.append(
            f"| {num} | `{acc}` | {desc} | {q_cov:.1f}% | {pct_id:.2f}% | {evalue:.2e} | {bit_score:.1f} |"
        )

    with open("annotation_table.md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print("Saved 'blast_output.xml' and 'annotation_table.md'.")


if __name__ == "__main__":
    q_header, q_seq = read_fasta(QUERY_FILE)
    rid, rtoe = submit_blast(q_seq)
    xml_data = fetch_results(rid, rtoe)
    parse_and_save(xml_data, q_header, len(q_seq))
