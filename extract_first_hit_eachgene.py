import os

name_file = "hq.name"#hq.name: MAG\ttaxonomy information
out_file = "merged_top_hit.tsv"

with open(name_file) as f:
   lines = [line.strip() for line in f if line.strip()]


with open(out_file, "w") as out:
    #out.write("query_with_sample\tsseqid\ttaxid\n")

    for line in lines:
        parts = line.split(",")

        faa_name = parts[0]

        # ✅ 提取目录名（第一个下划线前）
        sample_dir = faa_name.split('_')[0]

        # ✅ 提取你要加的前缀（去掉 .fa_bin.faa）
        prefix = faa_name.replace(".fa_bin.faa", "")

        diamond_file = os.path.join(sample_dir, "max300_subseq.fmt6")  # 按实际文件名改

        if not os.path.exists(diamond_file):
            print(f"Missing: {diamond_file}")
            continue

        seen = set()   # 每个样本独立

        print(f'processing:{faa_name}')
        with open(diamond_file) as f:
            for line in f:
                cols = line.strip().split("\t")

                qseqid = cols[0]
                sseqid = cols[1]
                taxid  = cols[-1]   # 最后一列

                if qseqid in seen:
                    continue
                seen.add(qseqid)

                new_query = f"{faa_name}_{qseqid[:-3]}"


                out.write(f"{new_query}\t{sseqid}\t{taxid}\n")
