#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import json
import os
import time
from collections import defaultdict

import requests

MAPPING_FILE = "idmapping_2025_12_02.tsv"
OUT_DIR = "cancer_distribution_tsv"
GENE_PREFIX = "ENSG"

# ----------------------------------------------------
# 1. 读取 Ensembl 基因 ID（To 列），去掉版本号、去重
# ----------------------------------------------------
def load_gene_ids(mapping_file):
    gene_ids = set()

    with open(mapping_file, "r", newline="") as f:
        # 先读表头，找到哪一列是 To（防止叫 "To " 之类）
        header_line = f.readline().rstrip("\n")
        headers = header_line.split("\t")
        # 找到名字中去掉空格后等于 "to" 的那一列
        try:
            to_idx = [i for i, h in enumerate(headers) if h.strip().lower() == "to"][0]
        except IndexError:
            raise RuntimeError("没有找到名为 'To' 的列，请检查表头")

        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) <= to_idx:
                continue
            raw_id = cols[to_idx].strip()
            if not raw_id:
                continue
            # 去掉版本号：ENSG00000100342.3 -> ENSG00000100342
            gene_id = raw_id.split(".")[0]
            # 只保留 Ensembl gene ID
            if gene_id.startswith(GENE_PREFIX):
                gene_ids.add(gene_id)

    return sorted(gene_ids)


# ----------------------------------------------------
# 2. 对单个基因：从 ssm_occurrences 拿到 project -> {case_ids}
# ----------------------------------------------------
SSM_OCC_ENDPOINT = "https://api.gdc.cancer.gov/ssm_occurrences"


def get_mutated_cases_by_project(gene_id, session):
    """
    返回: dict[project_id] = set(case_id)
    """
    fields = ["case.case_id", "case.project.project_id"]
    params = {
        "fields": ",".join(fields),
        "size": 100000,  # 足够大，一次性拿完
        "filters": json.dumps(
            {
                "op": "and",
                "content": [
                    {
                        "op": "in",
                        "content": {
                            # 这个字段名来自 GDC-QAG 的示例代码
                            # ssm.consequence.transcript.gene.symbol/gene_id
                            "field": "ssm.consequence.transcript.gene.gene_id",
                            "value": [gene_id],
                        },
                    }
                ],
            }
        ),
    }

    resp = session.get(SSM_OCC_ENDPOINT, params=params)
    if not resp.ok:
        raise RuntimeError(
            f"ssm_occurrences 请求失败 ({resp.status_code}): {resp.text[:500]}"
        )

    data = resp.json()["data"]["hits"]
    proj_to_cases = defaultdict(set)

    for item in data:
        proj = item["case"]["project"]["project_id"]
        case_id = item["case"]["case_id"]
        if proj and case_id:
            proj_to_cases[proj].add(case_id)

    return proj_to_cases


# ----------------------------------------------------
# 3. 对 project：用 case_ssms 算“有 SSM 数据的总 case 数”
# ----------------------------------------------------
CASE_SSMS_ENDPOINT = "https://api.gdc.cancer.gov/case_ssms"


def get_total_ssm_cases_for_project(project_id, session, cache):
    """
    返回指定 project 中 available_variation_data 包含 'ssm' 的 case 数。
    结果会缓存到 cache[project_id] 里。
    """
    if project_id in cache:
        return cache[project_id]

    params = {
        "fields": "project.project_id,available_variation_data",
        "size": 1,  # 只看 pagination.total，不需要真正的列表
        "filters": json.dumps(
            {
                "op": "and",
                "content": [
                    {
                        "op": "in",
                        "content": {
                            "field": "available_variation_data",
                            "value": ["ssm"],
                        },
                    },
                    {
                        "op": "=",
                        "content": {
                            "field": "project.project_id",
                            "value": project_id,
                        },
                    },
                ],
            }
        ),
    }

    resp = session.get(CASE_SSMS_ENDPOINT, params=params)
    if not resp.ok:
        raise RuntimeError(
            f"case_ssms 请求失败 ({resp.status_code}): {resp.text[:500]}"
        )

    total = resp.json()["data"]["pagination"]["total"]
    cache[project_id] = total
    return total


# ----------------------------------------------------
# 4. 写单个基因的 TSV
# ----------------------------------------------------
def write_cancer_distribution_tsv(gene_id, proj_to_cases, proj_total_ssm_cases):
    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = os.path.join(OUT_DIR, f"{gene_id}_cancer_distribution_ssm.tsv")

    rows = []
    for proj, case_set in proj_to_cases.items():
        cases_affected = len(case_set)
        total_ssm = proj_total_ssm_cases.get(proj, 0)
        pct = (cases_affected / total_ssm * 100) if total_ssm else None
        rows.append(
            {
                "project_id": proj,
                "cases_affected": cases_affected,
                "total_cases_with_ssm_data": total_ssm,
                "percent_cases_affected": pct,
            }
        )

    # 按百分比从大到小排序
    rows.sort(key=lambda r: (r["percent_cases_affected"] or 0), reverse=True)

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "project_id",
                "cases_affected",
                "total_cases_with_ssm_data",
                "percent_cases_affected",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    return out_path


# ----------------------------------------------------
# 5. 主程序：逐个基因执行，并实时打印进度
# ----------------------------------------------------
def main():
    print(f"读取 Ensembl 基因 ID 自 {MAPPING_FILE} ...")
    gene_ids = load_gene_ids(MAPPING_FILE)
    print(f"共找到 {len(gene_ids)} 个唯一 Ensembl 基因 ID")

    session = requests.Session()
    proj_total_cache = {}

    for idx, gene_id in enumerate(gene_ids, start=1):
        print(f"\n[{idx}/{len(gene_ids)}] 处理基因 {gene_id} ...", flush=True)

        try:
            proj_to_cases = get_mutated_cases_by_project(gene_id, session)
        except Exception as e:
            print(f"  ✗ 获取突变 case 失败：{e}")
            continue

        if not proj_to_cases:
            print("  ▹ 没有找到任何带该基因突变的 case，跳过（不生成 TSV）")
            continue

        total_occ = sum(len(s) for s in proj_to_cases.values())
        print(
            f"  ▹ 找到 {len(proj_to_cases)} 个 project，共 {total_occ} 个 "
            f"mutated case（去重前 occurrences）"
        )

        # 对每个 project 拿 denominator（有 SSM 数据的 case 总数）
        proj_total_ssm_cases = {}
        for proj in proj_to_cases.keys():
            try:
                total_ssm = get_total_ssm_cases_for_project(
                    proj, session, proj_total_cache
                )
            except Exception as e:
                print(f"    ✗ 获取 project {proj} 的总 SSM case 数失败：{e}")
                total_ssm = 0
            proj_total_ssm_cases[proj] = total_ssm
            # 稍微 sleep 一下，避免太频繁
            time.sleep(0.1)

        out_path = write_cancer_distribution_tsv(
            gene_id, proj_to_cases, proj_total_ssm_cases
        )
        print(f"  ✓ 已保存 TSV -> {out_path}", flush=True)

        # 根据需要可以再 sleep 一点，降低对 GDC 的压力
        time.sleep(0.2)


if __name__ == "__main__":
    main()
