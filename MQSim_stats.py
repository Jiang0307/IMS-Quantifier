import subprocess  
import re
import glob
import os

def run_mqsim_and_collect_stats():
    process = subprocess.Popen(
        [".\\MQSim.exe", "-i", ".\\SSD_config.xml", "-w", ".\\workload.xml"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        stdin=subprocess.PIPE,
        text=True,
        bufsize=1
    )

    total_write_us = 0
    total_read_us = 0
    deleted_files = set()

    for line in process.stdout:
        print(line, end="")

        # 自動送出 Enter
        if "Simulation complete; Press any key to exit." in line:
            process.stdin.write("\n")
            process.stdin.flush()

        # 擷取統計數據
        write_match = re.search(r"Total write time\s*:\s*(\d+)us", line)
        read_match = re.search(r"Total read time\s*:\s*(\d+)us", line)

        if write_match:
            total_write_us += int(write_match.group(1))
        if read_match:
            total_read_us += int(read_match.group(1))

        # 立即刪除新產生的 workload_scenario_*.xml
        for file in glob.glob("workload_scenario_*.xml"):
            if file not in deleted_files:
                try:
                    os.remove(file)
                    deleted_files.add(file)
                except Exception as e:
                    None

    process.wait()

    print("\n------------------ RESULT ------------------")
    print(f" - Total write time : {total_write_us} us = {total_write_us / 1_000_000:.4f} s")
    print(f" - Total read time  : {total_read_us} us = {total_read_us / 1_000_000:.4f} s")
    print(f" - Total runtime    : {total_read_us + total_write_us} us = {(total_read_us + total_write_us) / 1_000_000:.4f} s\n")

if __name__ == "__main__":
    run_mqsim_and_collect_stats()