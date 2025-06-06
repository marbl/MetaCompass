import sys
import re
from pathlib import Path

def get_reference_list(fpath):
    ref_names = []
    selected_ref_names = []
    with open(fpath, "r") as f:

        for l in f:
            ref_name_pre = l.split(".")[0]
            ref_name_post = l.split(".")[1][0]
            ref_name = f"{ref_name_pre}.{ref_name_post}"
            ref_names.append(ref_name)

            try:
                if re.sub(r'\s+', ' ', l).split(" ")[1] == "1":
                    selected_ref_names.append(ref_name)
            except Exception as e:
                print("Failed with error",e,"on line",re.sub(r'\s+', ' ', l).split(" "))

    return ref_names, selected_ref_names


def generate_trimmed_prokaryotes(ref_cov_path, new_file_name):
    prokaryotes_path = Path("/fs/cbcb-scratch/ayyangar/metacompass_new_version/tests/test_data/ref_db/prokaryotes.txt")
    new_file_path = Path(__file__).resolve().parent / new_file_name
    selected_new_file_path = Path(__file__).resolve().parent / f"selected_{new_file_name}"

    ref_names, selected_ref_names = get_reference_list(ref_cov_path)

    all_lines = []
    with open(prokaryotes_path, "r") as p:
        for l in p:
            all_lines.append(l.strip())

    content = [all_lines[0]]
    selected_content = [all_lines[0]]
    for ref_name in ref_names:
        for l in all_lines:
            if ref_name in l:
                content.append(l)
                break

    with open(new_file_path, "w") as trimmed_p:
        for line in content:
            trimmed_p.write(line + "\n")

    for ref_name in selected_ref_names:
        for l in all_lines:
            if ref_name in l:
                selected_content.append(l)
                break

    with open(selected_new_file_path, "w") as trimmed_p:
        for line in selected_content:
            trimmed_p.write(line + "\n")

ref_cov_path = Path(sys.argv[1])
new_file_name = sys.argv[2]

generate_trimmed_prokaryotes(ref_cov_path, new_file_name)