# 导入必要的库
from Bio.PDB import PDBParser, MMCIFParser, PDBList
import numpy as np
import pandas as pd
import os
import urllib.request
from tqdm import tqdm
from joblib import Parallel, delayed
from concurrent.futures import ThreadPoolExecutor, as_completed

# 定义 RNA/DNA 残基原子名称映射表（用于粗粒化坐标提取）
aa2longalt = {
    'DA': (" P  ", " O5'", " C5'", " C4'", " C3'", " O3'", " N9 ",
           " OP1", " OP2", " O4'", " C2'", " C1'",
           " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N6 ", None, None,
           "H5''", " H5'", " H4'", " H3'", "H2''", " H2'", " H1'", " H2 ", " H61", " H62", " H8 ", None),
    'DC': (" P  ", " O5'", " C5'", " C4'", " C3'", " O3'", " N1 ",
           " OP1", " OP2", " O4'", " C2'", " C1'",
           " C2 ", " O2 ", " N3 ", " C4 ", " N4 ", " C5 ", " C6 ", None, None, None, None,
           "H5''", " H5'", " H4'", " H3'", "H2''", " H2'", " H1'", " H42", " H41", " H5 ", " H6 ", None),
    'DG': (" P  ", " O5'", " C5'", " C4'", " C3'", " O3'", " N9 ",
           " OP1", " OP2", " O4'", " C2'", " C1'",
           " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N2 ", " O6 ", None,
           "H5''", " H5'", " H4'", " H3'", "H2''", " H2'", " H1'", " H1 ", " H22", " H21", " H8 ", None),
    'DT': (" P  ", " O5'", " C5'", " C4'", " C3'", " O3'", " N1 ",
           " OP1", " OP2", " O4'", " C2'", " C1'",
           " C2 ", " O2 ", " N3 ", " C4 ", " O4 ", " C5 ", " C7 ", " C6 ", None, None, None,
           "H5''", " H5'", " H4'", " H3'", "H2''", " H2'", " H1'", " H3 ", " H71", " H72", " H73", " H6 "),
    'DX': (" P  ", " O5'", " C5'", " C4'", " C3'", " O3'", None,
           " OP1", " OP2", " O4'", " C2'", " C1'",
           None, None, None, None, None, None, None, None, None, None, None,
           "H5''", " H5'", " H4'", " H3'", "H2''", " H2'", " H1'", None, None, None, None, None),
    'A': (" P  ", " O5'", " C5'", " C4'", " C3'", " O3'", " N9 ",
          " OP1", " OP2", " O4'", " C1'", " C2'", " O2'",
          " N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 ", " N6 ", " N7 ", " C8 ", None,
          " H5'", "H5''", " H4'", " H3'", " H2'", "HO2'", " H1'", " H2 ", " H61", " H62", " H8 ", None),
    'C': (" P  ", " O5'", " C5'", " C4'", " C3'", " O3'", " N1 ",
          " OP1", " OP2", " O4'", " C1'", " C2'", " O2'",
          " C2 ", " O2 ", " N3 ", " C4 ", " N4 ", " C5 ", " C6 ", None, None, None,
          " H5'", "H5''", " H4'", " H3'", " H2'", "HO2'", " H1'", " H42", " H41", " H5 ", " H6 ", None),
    'G': (" P  ", " O5'", " C5'", " C4'", " C3'", " O3'", " N9 ",
          " OP1", " OP2", " O4'", " C1'", " C2'", " O2'",
          " N1 ", " C2 ", " N2 ", " N3 ", " C4 ", " C5 ", " C6 ", " O6 ", " N7 ", " C8 ",
          " H5'", "H5''", " H4'", " H3'", " H2'", "HO2'", " H1'", " H1 ", " H22", " H21", " H8 ", None),
    'U': (" P  ", " O5'", " C5'", " C4'", " C3'", " O3'", " N1 ",
          " OP1", " OP2", " O4'", " C1'", " C2'", " O2'",
          " C2 ", " O2 ", " N3 ", " C4 ", " O4 ", " C5 ", " C6 ", None, None, None,
          " H5'", "H5''", " H4'", " H3'", " H2'", "HO2'", " H1'", " H3 ", " H5 ", " H6 ", None, None),
    'RX': (" P  ", " O5'", " C5'", " C4'", " C3'", " O3'", None,
           " OP1", " OP2", " O4'", " C1'", " C2'", " O2'",
           None, None, None, None, None, None, None, None, None, None,
           " H5'", "H5''", " H4'", " H3'", " H2'", "HO2'", " H1'", None, None, None, None, None)
}

# 打印支持的残基类型（调试用）
print(aa2longalt.keys())

# 下载 PDB 文件（优先 .pdb，失败则尝试 .cif）
def download_pdb(pdb_id, pdb_out_dir):
    pdb_path = os.path.join(pdb_out_dir, f"{pdb_id}.pdb")
    cif_path = os.path.join(pdb_out_dir, f"{pdb_id}.cif")

    if os.path.exists(pdb_path) or os.path.exists(cif_path):
        return True

    try:
        url = f"http://files.rcsb.org/download/{pdb_id}.pdb"
        urllib.request.urlretrieve(url, pdb_path)
        return True
    except Exception as e1:
        try:
            mmcif_url = f"http://files.rcsb.org/download/{pdb_id}.cif"
            urllib.request.urlretrieve(mmcif_url, cif_path)
            return True
        except Exception as e2:
            print(f"Error downloading {pdb_id}: {e1}, {e2}")
            return False

# 下载对应的 PDB 文件
def load_and_filter(pdb_list, pdb_out_dir, n_jobs=4):
    if not os.path.exists(pdb_out_dir):
        os.makedirs(pdb_out_dir)

    print(f"Processing {len(pdb_list)} items")


    # 获取唯一 PDB ID 列表
    pdb_ids = list(set(pdb_list))

    # 多线程下载
    results = Parallel(n_jobs=n_jobs)(
        delayed(download_pdb)(pdb_id, pdb_out_dir) for pdb_id in tqdm(pdb_ids)
    )

    success_count = sum(results)
    print(f"Successfully saved {success_count} out of {len(pdb_ids)} complexes in {pdb_out_dir}")

# 提取 RNA 序列和粗粒化坐标（C3' 坐标）
def get_rna_seq_and_coords(pdb_id, pdb_dir='raw_pdb'):
    try:
        # 尝试加载本地文件
        if os.path.exists(os.path.join(pdb_dir, f"{pdb_id}.pdb")):
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(pdb_id, os.path.join(pdb_dir, f"{pdb_id}.pdb"))
        elif os.path.exists(os.path.join(pdb_dir, f"{pdb_id}.cif")):
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure(pdb_id, os.path.join(pdb_dir, f"{pdb_id}.cif"))
        else:
            raise FileNotFoundError(f"{pdb_id} not found locally.")
    except Exception as e:
        try:
            # 若本地不存在，尝试从 RCSB 下载
            pdb_list = PDBList()
            pdb_list.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=pdb_dir)
            if os.path.exists(os.path.join(pdb_dir, f"{pdb_id}.pdb")):
                parser = PDBParser(PERMISSIVE=1, QUIET=True)
                structure = parser.get_structure(pdb_id, os.path.join(pdb_dir, f"{pdb_id}.pdb"))
            else:
                raise FileNotFoundError(f"{pdb_id} not found after download.")
        except Exception as e2:
            print(f"Error processing {pdb_id}: {e}, {e2}")
            return None, None, None

    seq = ""
    coords = []
    chain_ids = []

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = str(residue.resname)
                if resname in aa2longalt:
                    res_coord = []
                    for atomname in aa2longalt[resname]:
                        if atomname is None:
                            res_coord.append(np.array([np.nan] * 3))
                        elif atomname.replace(" ", "") not in residue:
                            res_coord.append(np.array([np.nan] * 3))
                        else:
                            res_coord.append(residue[atomname.replace(" ", "")].coord)
                    coords.append(res_coord)
                    seq += resname[-1]  # 取最后一个字符作为单字母表示
                    chain_ids.append(chain.id)
        break  # 只处理第一个 model

    assert len(seq) == len(coords), f"Length mismatch for {pdb_id}"

    if seq == "":
        return None, None, None

    return seq, np.array(coords), chain_ids

# 处理单个 PDB 文件：保存序列和坐标
def process_single_pdb(pdb_id):
    seq, coords, chain_ids = get_rna_seq_and_coords(pdb_id)
    if seq is None:
        return

    chain_id_set = set(chain_ids)
    cid_array = np.array(chain_ids)
    seq_arr = np.array(list(seq))

    for cid in chain_id_set:
        mask = (cid_array == cid)
        chain_seq = ''.join(seq_arr[mask])
        chain_coords = coords[mask]
        chain_c3_coords = chain_coords[:, 4, :]  # 第5个原子是 C3'

        os.makedirs('./seq', exist_ok=True)
        os.makedirs('./C3_coords', exist_ok=True)
        os.makedirs('./coords', exist_ok=True)

        with open(f'./seq/{pdb_id}_{cid}.fasta', 'w') as f:
            f.write(f'>{pdb_id}_{cid}\n{chain_seq}\n')

        np.save(f'./C3_coords/{pdb_id}_{cid}.npy', chain_c3_coords)
        np.save(f'./coords/{pdb_id}_{cid}.npy', chain_coords)
