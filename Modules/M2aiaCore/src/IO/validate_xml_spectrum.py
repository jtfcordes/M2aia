import argparse
import xml.etree.ElementTree as ET
from pathlib import Path

NAMESPACE = {"ns": "http://psi.hupo.org/ms/mzml"}

# Prepend namespace to the element
def ns_tag(tag):
    return f"{{{NAMESPACE['ns']}}}{tag}"

def parse_arguments():
    parser = argparse.ArgumentParser(description="Validate and fix XML imzML file structure.")
    parser.add_argument("input_file", type=str, help="Path to the XML file to be validated.")
    parser.add_argument("--dryrun", action="store_true", help="Enable dry-run mode to list findings and possible corrections.")
    parser.add_argument("--output_file", type=str, help="Path to save the corrected XML file (if dry-run is disabled).")
    return parser.parse_args()

def check_and_correct_file_content(root, dryrun):
    # Define the required cvParam elements
    ims_params ={
        "IMS:1000031": "processed",
        "IMS:1000030": "continuous"
    }

    ms_params = {
        "MS:1000127": "centroid spectrum",
        "MS:1000128": "profile spectrum"
    }

    # Scan the entire XML for IMS and MS parameters
    all_ims_accessions = set()
    all_ms_accessions = set()
    for elem in root.findall(".//ns:cvParam", NAMESPACE):
        if elem.attrib.get("cvRef") == "IMS":
            all_ims_accessions.add(elem.attrib.get("accession"))
        elif elem.attrib.get("cvRef") == "MS":
            all_ms_accessions.add(elem.attrib.get("accession"))

    # Find <fileContent> section
    file_content = root.find(".//ns:fileContent", NAMESPACE)
    if file_content is None:
        print("[ERROR] <fileContent> tag is missing!")
        return False

    # Check for existing cvParam elements in <fileContent>
    existing_ims = set()
    existing_ms = set()
    for elem in file_content.findall("ns:cvParam", NAMESPACE):
        if elem.attrib.get("cvRef") == "IMS":
            existing_ims.add(elem.attrib.get("accession"))
        elif elem.attrib.get("cvRef") == "MS":
            existing_ms.add(elem.attrib.get("accession"))

    # Determine missing IMS and MS params
    # missing_ims = [p for p in ims_params if p["accession"] not in existing_ims]
    # missing_ms = [p for p in ms_params if p["accession"] not in existing_ms]

    # Report findings in dry-run mode
    if dryrun:
        print("All IMS parameters found:")
        print(f"  {all_ims_accessions}")
        print("All MS parameters found:")
        print(f"  {all_ms_accessions}")

        print("[DRY-RUN] Findings:")
        print(f"  IMS parameters found: {existing_ims}")
        found = [acc + ": " + "..." if not acc in ims_params else acc + ": " + ims_params[acc] for acc in existing_ims ]
        for param in found:
            print(f"  - {param}")
        print(f"  MS parameters found: {existing_ms}")
        print("[DRY-RUN] Missing in <fileContent>:")
        # for param in missing_ims + missing_ms:
        #     print(f"  - {param}")
        return True

    # Add missing parameters to <fileContent> if not dry-run
    for param in missing_ims + missing_ms:
        new_elem = ET.SubElement(file_content, ns_tag("cvParam"), param)
        print(f"[FIX] Added {param} to <fileContent>.")

    return True

def main():
    args = parse_arguments()
    input_file = Path(args.input_file)

    if not input_file.is_file():
        print(f"[ERROR] File not found: {input_file}")
        return

    # Parse the XML file
    try:
        tree = ET.parse(input_file)
        root = tree.getroot()
    except ET.ParseError as e:
        print(f"[ERROR] Failed to parse XML: {e}")
        return

    # Validate and optionally fix the XML
    success = check_and_correct_file_content(root, args.dryrun)
    if not success:
        return

    # Save the corrected XML if not dry-run
    if not args.dryrun:
        output_file = Path(args.output_file or "corrected_" + input_file.name)
        tree.write(output_file, encoding="utf-8", xml_declaration=True)
        print(f"[INFO] Corrected XML saved to {output_file}")

if __name__ == "__main__":
    main()
