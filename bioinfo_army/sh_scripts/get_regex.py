import argparse as ap
import xml.etree.ElementTree as ET

def read_xml(xml_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    return root

def find_regex(root):
    motifs = root.findall(".//motif")
    regex = motifs[0].findtext("regular_expression").strip()
    return regex

def main():
    parser = ap.ArgumentParser(description="Reads a meme.xml file and outputs the regular expression (regex)")
    parser.add_argument("file", help="Path of the meme.xml file",type=str)
    args = parser.parse_args()
    root = read_xml(args.file)
    regex = find_regex(root)
    print(regex, end="")

if __name__ == "__main__":
    main()

