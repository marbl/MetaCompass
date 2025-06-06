from typing import List
from pathlib import Path
import subprocess


def run_shell_cmd(cmds: List, outfile=None):
    """
    Run one or more shell commands, optionally piping them together.

    Args
        cmds (List): A single command as a list or multiple commands as a list of lists.
        outfile (Path): Optional output file to redirect the final command's output.
    """
    # print("Command", cmds)

    try:
        if not isinstance(cmds[0], list):  # If cmds is a single command
            cmds = [cmds]

        # Start the first process
        processes = [subprocess.Popen(cmds[0], stdout=subprocess.PIPE)]

        # Start the rest of the processes, connecting each one's input to the previous one's output
        for cmd in cmds[1:]:
            processes.append(subprocess.Popen(cmd, stdin=processes[-1].stdout, stdout=subprocess.PIPE))
            processes[-2].stdout.close()  # Allow previous process to receive a SIGPIPE if the next process exits

        # Wait for the final process to complete and capture its output
        final_process = processes[-1]
        stdout, stderr = final_process.communicate()

        # Check for errors
        if final_process.returncode != 0:
            raise subprocess.CalledProcessError(final_process.returncode, final_process.args)

        # Write the output to a file if given
        if outfile:
            with open(outfile, 'a') as f:
                f.write(stdout.decode())

        return stdout.decode()

    except subprocess.CalledProcessError as e:
        error_message = e.stderr.decode() if e.stderr else "No output"
        print(f"Command {e.cmd} failed with return code {e.returncode} - {error_message}")
    except Exception as e:
        print(f"An error occurred - {str(e)}.")


def run_metacompass(metacompass_runner, metacompass_path, reference_db, result_output_path,
                    forward_read_path, reverse_read_path):
    cmd = [metacompass_runner.resolve().as_posix(),
           result_output_path, metacompass_path, reference_db, forward_read_path, reverse_read_path]

    subprocess.run(cmd)


def read_file(fpath: Path):
    lines = []
    with open(fpath, "r") as f:
        lines = [line.strip() for line in f.readlines()]

    return lines


def test_reference_selection_and_culling(tmp_dir, sel_cull_ref_db, metacompass_path, metacompass_runner):
    ref_db = sel_cull_ref_db.ref_db
    forward_read = sel_cull_ref_db.forward_read
    reverse_read = sel_cull_ref_db.reverse_read
    result_output_path = tmp_dir / "reference_db_test"

    run_metacompass(
        metacompass_runner=metacompass_runner,
        metacompass_path=metacompass_path.resolve().as_posix(),
        reference_db=ref_db.resolve().as_posix(),
        result_output_path=result_output_path.resolve().as_posix(),
        forward_read_path=forward_read.resolve().as_posix(),
        reverse_read_path=reverse_read.resolve().as_posix()
    )

    actual_ref_cul = result_output_path / "reference_culling"
    actual_ref_sel = result_output_path / "reference_selection"

    # check if GCA_009556455.1_ASM955645v1_genomic.fna was selected

    ref_selected = read_file(actual_ref_sel / "cluster_refs" / "reference_candidates.txt")
    assert len(ref_selected) == 1
    assert ref_selected[0] == "GCA_009556455.1"

    # tests for culling

    cluster_content = read_file(actual_ref_cul / "cluster_content.csv")
    assert len(cluster_content) == 1  # only 1 cluster
    assert cluster_content[0].split(",")[1].split("/")[-1] == "GCA_009556455.1_ASM955645v1_genomic.fna"

    cluster_1_fna = read_file(actual_ref_cul / "cluster_1.fna")
    assert len(cluster_1_fna) == 17040
