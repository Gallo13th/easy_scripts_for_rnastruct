import re
import os

RFAM_ROOT_DIR = "/Users/letiangao/rfam"

def extract_all_gf_fields(seed_file): 
    """
    提取所有 #=GF 字段的内容
    返回字典格式，键为字段名，值为列表（因为同一字段可能有多个值）
    """
    extracted = {}
    
    lines = open(seed_file,'r').readlines()
    for line in lines:
        # 匹配 #=GF 开头的行
        if line.startswith('#=GF'):
            # 使用正则表达式分割，最多分割成3部分：#=GF、字段名、值
            parts = re.split(r'\s+', line, maxsplit=2)
            if len(parts) >= 3:
                tag = parts[1]      # 字段名（如 AC, ID, DE 等）
                value = parts[2]    # 字段值
                
                # 如果该字段名还未在字典中，初始化为空列表
                if tag not in extracted:
                    extracted[tag] = []
                
                # 将值添加到对应字段的列表中
                extracted[tag].append(value)
    for key, value in extracted.items():
        extracted[key] = ''.join(value).rstrip().replace('\n',';')    
    return extracted

def print_gf_fields(gf_dict):
    """
    格式化打印提取的 GF 字段
    """
    for tag, values in gf_dict.items():
        print(f"#=GF {tag}:")
        for value in values:
            print(f"  {value}")
        print()  # 空行分隔

def main():
    import pandas as pd
    from tqdm import tqdm
    results = []
    for rfam_dir in tqdm(os.listdir(RFAM_ROOT_DIR)):
        if not rfam_dir.startswith("RF"):
            continue
        seed_path = os.path.join(RFAM_ROOT_DIR,rfam_dir,"Rfam.seed")
        results.append(extract_all_gf_fields(seed_path))
    df = pd.DataFrame(results)
    print(df)
    df.to_csv('./rfam/rfam_seed_header.csv',sep='\t',index=False)
    
        
if __name__ == '__main__':
    main()