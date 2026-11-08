#!/usr/bin/env python3
"""
RSV Whole Genome Sequence Downloader from GenBank
==================================================

This script downloads all available RSV (Respiratory Syncytial Virus) 
whole genome sequences from NCBI GenBank using multiple approaches.

阅读这份代码时，先抓住主流程：
1. search_rsv_genomes：搜索候选 accession number（登录号）
2. get_sequence_info：取得每个 accession 的摘要 metadata
3. filter_rsv_genomes：筛选看起来像 RSV 完整基因组的记录
4. download_sequences：下载 FASTA 和 GenBank 两种格式
5. save_metadata：保存 CSV、JSON 和文字摘要
6. create_combined_fasta：把单独的 FASTA 合并成一个文件
7. analyze_dataset：计算数量、长度、年份等简单统计

注意：这段程序会访问 NCBI，并可能创建上万个文件。阅读和理解后再运行。

Requirements:
- biopython (pip install biopython)
- pandas (optional, for data analysis)
- requests (for direct NCBI API access)

Author: Claude Assistant
Date: 2026
"""

# =========================
# 1. 导入程序所需的工具
# =========================
# os：创建文件夹、拼接路径、检查文件是否已经存在。
import os
# sys：把日志输出到终端，并在缺少依赖时退出程序。
import sys
# time：在 NCBI 请求之间暂停，避免请求过快。
import time
# json / csv：把 metadata 保存为两种常用表格/数据格式。
import json
import csv
# datetime：把实际下载时间写入 summary 文件。
from datetime import datetime
# 这些类型标记帮助读者理解参数和返回值；不会改变程序运行方式。
# 本版本中 Optional 和 Tuple 实际上没有被使用。
from typing import List, Dict, Optional, Tuple
# logging：同时在终端和日志文件中记录进度及错误。
import logging

# =========================
# 2. 设置运行日志 logging
# =========================
# INFO 级别会显示主要进度；DEBUG 信息默认不会显示。
# 每条日志同时写入：
#   ① 当前工作目录下的 rsv_download.log
#   ② 终端屏幕（stdout）
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('rsv_download.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# =========================
# 3. 检查可用的下载工具
# =========================
# 首选 Biopython，因为 Entrez 模块专门封装了 NCBI E-utilities。
# 如果没有安装 Biopython，程序会尝试用 requests 直接调用 NCBI API。
try:
    from Bio import Entrez, SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    # 注：SeqIO、Seq、SeqRecord 在本版本中没有实际使用；真正使用的是 Entrez。
    BIOPYTHON_AVAILABLE = True
except ImportError:
    logger.warning("Biopython not available. Some features will be limited.")
    BIOPYTHON_AVAILABLE = False

# requests 是备用网络访问方式：Biopython 不存在时才会使用。
try:
    import requests
    REQUESTS_AVAILABLE = True
except ImportError:
    logger.warning("Requests library not available. Direct API access will be limited.")
    REQUESTS_AVAILABLE = False

# pandas 被标记为“可选”，但本版本后面的代码实际上没有使用 pd。
# 因此缺少 pandas 不会阻止下载。
try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    logger.warning("Pandas not available. Data analysis features will be limited.")
    PANDAS_AVAILABLE = False


# ================================================
# 4. 下载器类：把共享的设置和各个步骤组织在一起
# ================================================
class RSVGenBankDownloader:
    """负责搜索、筛选、下载并汇总 RSV GenBank 记录的对象。"""
    
    def __init__(self, email: str, output_dir: str = "rsv_genomes"):
        """
        Initialize the downloader
        
        Args:
            email: Your email address (required by NCBI)
            output_dir: Directory to save downloaded sequences
        """
        # self 保存“整次运行都要使用”的信息，后面的每个 method 都能访问。
        self.email = email
        self.output_dir = output_dir
        # NCBI E-utilities 的共同网址；后面会加上 esearch、esummary 或 efetch。
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        
        # exist_ok=True：文件夹已存在时不报错，因此程序可以继续未完成的下载。
        os.makedirs(output_dir, exist_ok=True)
        
        if BIOPYTHON_AVAILABLE:
            # NCBI 要求自动请求提供联系邮箱；若发生滥用或故障可联系使用者。
            Entrez.email = email
            
        logger.info(f"Initialized RSV GenBank Downloader")
        logger.info(f"Email: {email}")
        logger.info(f"Output directory: {output_dir}")
    
    def search_rsv_genomes(self, retmax: int = 10000) -> List[str]:
        """
        Search for RSV whole genome sequences in GenBank
        
        Args:
            retmax: Maximum number of results to return
            
        Returns:
            List of GenBank accession IDs
        """
        logger.info("Searching for RSV whole genome sequences...")
        
        # 一个搜索式可能漏掉命名方式不同的记录，所以使用四个有重叠的搜索式。
        # [Organism]、[All Fields]、[SLEN]、[Title] 是 NCBI 搜索字段标签。
        # 15000:16000[SLEN] 表示序列长度在约 15–16 kb 之间。
        search_terms = [
            '("respiratory syncytial virus"[Organism] AND complete genome)',
            '("human respiratory syncytial virus"[Organism] AND complete genome)',
            '("RSV"[All Fields] AND complete genome AND 15000:16000[SLEN])',  # RSV genome size
            '("respiratory syncytial virus"[Organism] AND "complete genome"[Title])',
        ]
        
        # set 会自动去重。同一 accession 可能被多个搜索式找到，但只保留一次。
        all_ids = set()
        
        # 依次运行每个搜索式，并把结果加入同一个集合。
        for search_term in search_terms:
            try:
                if BIOPYTHON_AVAILABLE:
                    # 首选路径：Biopython Entrez。
                    ids = self._search_with_biopython(search_term, retmax)
                elif REQUESTS_AVAILABLE:
                    # 备用路径：直接发送 HTTP 请求。
                    ids = self._search_with_requests(search_term, retmax)
                else:
                    logger.error("Neither Biopython nor requests available for search")
                    return []
                
                # update 会加入新 accession，并自动忽略已经存在的 accession。
                all_ids.update(ids)
                logger.info(f"Found {len(ids)} sequences with search: {search_term[:50]}...")
                # 主动暂停，降低对 NCBI 公共服务器的压力。
                time.sleep(0.5)
                
            except Exception as e:
                # 一个搜索式失败时记录错误，但继续尝试其余搜索式。
                logger.error(f"Error in search '{search_term}': {e}")
                continue
        
        # set 本身没有固定顺序，所以这里产生的 accession 顺序不保证每次相同。
        unique_ids = list(all_ids)
        logger.info(f"Total unique RSV genome sequences found: {len(unique_ids)}")
        return unique_ids
    
    def _search_with_biopython(self, search_term: str, retmax: int) -> List[str]:
        """使用 Biopython 调用 ESearch；只返回 accession，不下载序列。"""
        # db="nucleotide" 指定 NCBI Nucleotide 数据库。
        # retmax 是这一个搜索式最多返回多少条结果。
        # idtype="acc" 要求返回 accession，而不是内部数字 UID。
        handle = Entrez.esearch(
            db="nucleotide",
            term=search_term,
            retmax=retmax,
            idtype="acc"
        )
        # Entrez.read 把 NCBI 返回的数据解析成 Python 对象。
        search_results = Entrez.read(handle)
        # 网络 handle 使用后应关闭。
        handle.close()
        return search_results["IdList"]
    
    def _search_with_requests(self, search_term: str, retmax: int) -> List[str]:
        """不用 Biopython，直接请求 ESearch 的 JSON 结果。"""
        url = f"{self.base_url}esearch.fcgi"
        params = {
            'db': 'nucleotide',
            'term': search_term,
            'retmax': retmax,
            'retmode': 'json',
            'email': self.email
        }
        
        # params 会由 requests 安全地编码进 URL 查询参数。
        response = requests.get(url, params=params)
        # HTTP 4xx/5xx 时抛出异常，交给外层 try/except 记录。
        response.raise_for_status()
        data = response.json()
        return data.get('esearchresult', {}).get('idlist', [])
    
    def get_sequence_info(self, ids: List[str]) -> List[Dict]:
        """
        Get detailed information about sequences
        
        Args:
            ids: List of GenBank IDs
            
        Returns:
            List of dictionaries with sequence information
        """
        logger.info(f"Getting detailed information for {len(ids)} sequences...")
        
        # ESummary 返回“摘要 metadata”，不是核苷酸序列本身。
        # 每 100 个 accession 为一批，避免一次请求过大。
        batch_size = 100
        all_info = []
        
        # range(0, len(ids), 100) 会产生 0、100、200……这些批次起点。
        for i in range(0, len(ids), batch_size):
            # Python 切片取得当前批次；最后一批可以少于 100 条。
            batch_ids = ids[i:i+batch_size]
            logger.info(f"Processing batch {i//batch_size + 1}/{(len(ids)-1)//batch_size + 1}")
            
            try:
                if BIOPYTHON_AVAILABLE:
                    batch_info = self._get_info_with_biopython(batch_ids)
                elif REQUESTS_AVAILABLE:
                    batch_info = self._get_info_with_requests(batch_ids)
                else:
                    logger.error("No method available to get sequence info")
                    return []
                
                # extend 把这一批的多个 dictionary 加入总列表。
                all_info.extend(batch_info)
                time.sleep(1)  # 批次间暂停，尊重 NCBI 服务。
                
            except Exception as e:
                # 注意：整批失败时会跳过这一批；本代码没有自动重试。
                logger.error(f"Error processing batch {i//batch_size + 1}: {e}")
                continue
        
        return all_info
    
    def _get_info_with_biopython(self, ids: List[str]) -> List[Dict]:
        """用 Biopython ESummary 把 accession 转换成摘要 metadata。"""
        handle = Entrez.esummary(db="nucleotide", id=",".join(ids))
        summaries = Entrez.read(handle)
        handle.close()
        
        # 把 NCBI 字段整理成后续代码统一使用的较简单字段名。
        info_list = []
        for summary in summaries:
            info_list.append({
                'id': summary.get('AccessionVersion', ''),
                'title': summary.get('Title', ''),
                'length': summary.get('Length', 0),
                'organism': summary.get('Organism', ''),
                # create_date 是 GenBank 记录建立日期，不是样本采集日期。
                'create_date': summary.get('CreateDate', ''),
                'update_date': summary.get('UpdateDate', ''),
                'authors': summary.get('Authors', ''),
                'journal': summary.get('Source', '')
            })
        # 注意：ESummary 通常不提供 geo_loc_name 和 collection_date，
        # 所以这份 downloader 生成的 metadata CSV 没有国家/采样日期字段。
        return info_list
    
    def _get_info_with_requests(self, ids: List[str]) -> List[Dict]:
        """直接请求 ESummary JSON，并整理成与 Biopython 路径相同的结构。"""
        url = f"{self.base_url}esummary.fcgi"
        params = {
            'db': 'nucleotide',
            'id': ','.join(ids),
            'retmode': 'json',
            'email': self.email
        }
        
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()
        
        info_list = []
        # .get(..., default) 可在字段缺失时避免 KeyError。
        result = data.get('result', {})
        for uid in result.get('uids', []):
            if uid in result:
                summary = result[uid]
                info_list.append({
                    'id': summary.get('accessionversion', ''),
                    'title': summary.get('title', ''),
                    'length': summary.get('slen', 0),
                    'organism': summary.get('organism', ''),
                    'create_date': summary.get('createdate', ''),
                    'update_date': summary.get('updatedate', ''),
                    'authors': summary.get('authors', ''),
                    'journal': summary.get('source', '')
                })
        return info_list
    
    def filter_rsv_genomes(self, seq_info: List[Dict]) -> List[Dict]:
        """
        Filter sequences to only include RSV complete genomes
        
        Args:
            seq_info: List of sequence information dictionaries
            
        Returns:
            Filtered list of RSV genome sequences
        """
        logger.info("Filtering sequences for RSV complete genomes...")
        
        # 只有同时通过长度、RSV 名称、完整基因组描述和排除词检查的记录才保留。
        filtered = []
        for info in seq_info:
            # 转成小写后比较，可忽略原文字母大小写差异。
            title_lower = info['title'].lower()
            organism_lower = info['organism'].lower()
            length = int(info['length'])
            
            # 使用较宽的 14,500–16,000 nt 范围，容许天然长度差异和注释差异。
            is_right_length = 14500 <= length <= 16000
            
            # organism 或 title 只要包含任一 RSV 相关词，就通过名称检查。
            is_rsv = any(term in organism_lower or term in title_lower for term in [
                'respiratory syncytial virus',
                'rsv',
                'human respiratory syncytial'
            ])
            
            # 标题必须明确表示 complete/whole genome。
            is_complete = any(term in title_lower for term in [
                'complete genome',
                'complete sequence',
                'whole genome'
            ])
            
            # 标题含有这些词时排除，以减少 partial、单基因、载体等非整基因组记录。
            # 这是启发式规则（heuristic），不是绝对正确的生物学判定；
            # 例如某些完整基因组标题也可能出现 "gene"，从而被误排除。
            is_not_excluded = not any(term in title_lower for term in [
                'partial',
                'gene',
                'protein',
                'mrna',
                'clone',
                'vector'
            ])
            
            # 四个布尔条件必须全部为 True。
            if is_rsv and is_right_length and is_complete and is_not_excluded:
                filtered.append(info)
                logger.debug(f"Kept: {info['id']} - {info['title'][:60]}...")
            else:
                logger.debug(f"Filtered out: {info['id']} - {info['title'][:60]}...")
        
        logger.info(f"After filtering: {len(filtered)} RSV complete genomes")
        return filtered
    
    def download_sequences(self, seq_info: List[Dict]) -> None:
        """
        Download the actual sequence data
        
        Args:
            seq_info: List of sequence information dictionaries
        """
        logger.info(f"Downloading {len(seq_info)} RSV genome sequences...")
        
        # FASTA 主要保存序列；GenBank flat file 还保存 feature、country、date 等注释。
        # 两种格式分别放进 fasta/ 和 genbank/ 子文件夹。
        fasta_dir = os.path.join(self.output_dir, "fasta")
        genbank_dir = os.path.join(self.output_dir, "genbank")
        os.makedirs(fasta_dir, exist_ok=True)
        os.makedirs(genbank_dir, exist_ok=True)
        
        # 计数器用于最后的日志，不直接决定 metadata 中的记录数。
        successful_downloads = 0
        failed_downloads = 0
        
        # enumerate 同时提供从 0 开始的索引 i 和当前 metadata dictionary。
        for i, info in enumerate(seq_info):
            acc_id = info['id']
            logger.info(f"Downloading {i+1}/{len(seq_info)}: {acc_id}")
            
            try:
                # 同一 accession 下载两次：一次 FASTA，一次 GenBank。
                # 只有两次都没有抛出异常，successful_downloads 才加 1。
                self._download_single_sequence(acc_id, fasta_dir, "fasta")
                self._download_single_sequence(acc_id, genbank_dir, "gb")
                successful_downloads += 1
                
            except Exception as e:
                # 某个 accession 失败不会终止整个任务；程序继续下一个。
                # 但本函数不自动重试，也不会删除可能已成功写入的一种格式。
                logger.error(f"Failed to download {acc_id}: {e}")
                failed_downloads += 1
                continue
            
            # 每个 accession 后暂停。上万条记录会因此花费很长时间。
            time.sleep(0.5)
        
        logger.info(f"Download complete: {successful_downloads} successful, {failed_downloads} failed")
    
    def _download_single_sequence(self, acc_id: str, output_dir: str, format_type: str) -> None:
        """下载一个 accession 的一种格式，并写入一个独立文件。"""
        # 例如 PP376595.1 + fasta -> PP376595.1.fasta。
        filename = f"{acc_id}.{format_type}"
        filepath = os.path.join(output_dir, filename)
        
        # 已存在就跳过，允许中断后继续；但这里只检查“存在”，不验证文件是否完整。
        if os.path.exists(filepath):
            logger.debug(f"Skipping {acc_id} - already exists")
            return
        
        if BIOPYTHON_AVAILABLE:
            # EFetch 才真正返回序列/GenBank flat file 内容。
            handle = Entrez.efetch(
                db="nucleotide",
                id=acc_id,
                rettype=format_type,
                retmode="text"
            )
            sequence_data = handle.read()
            handle.close()
        elif REQUESTS_AVAILABLE:
            # 备用路径请求同一个 EFetch endpoint。
            url = f"{self.base_url}efetch.fcgi"
            params = {
                'db': 'nucleotide',
                'id': acc_id,
                'rettype': format_type,
                'retmode': 'text',
                'email': self.email
            }
            response = requests.get(url, params=params)
            response.raise_for_status()
            sequence_data = response.text
        else:
            raise Exception("No method available to download sequences")
        
        # 'w' 会创建新文件；若文件已存在则覆盖。
        # 上面的 exists 检查使正常流程不会覆盖已有文件。
        with open(filepath, 'w') as f:
            f.write(sequence_data)
        
        logger.debug(f"Saved {filepath}")
    
    def save_metadata(self, seq_info: List[Dict]) -> None:
        """把筛选后的摘要 metadata 保存成 JSON、CSV 和文字 summary。"""
        logger.info("Saving metadata...")
        
        # JSON 保留 list + dictionary 结构，适合程序读取。
        json_path = os.path.join(self.output_dir, "rsv_metadata.json")
        with open(json_path, 'w') as f:
            json.dump(seq_info, f, indent=2)
        
        # CSV 是平面表格，适合 Excel、Google Sheets 或 pandas。
        csv_path = os.path.join(self.output_dir, "rsv_metadata.csv")
        if seq_info:
            # 第一条记录的 keys 决定 CSV 列名及顺序。
            fieldnames = seq_info[0].keys()
            with open(csv_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(seq_info)
        
        # summary 是给人快速查看的文字文件。
        summary_path = os.path.join(self.output_dir, "download_summary.txt")
        with open(summary_path, 'w') as f:
            f.write(f"RSV GenBank Download Summary\n")
            f.write(f"============================\n\n")
            f.write(f"Download Date: {datetime.now()}\n")
            f.write(f"Total Sequences: {len(seq_info)}\n")
            f.write(f"Output Directory: {self.output_dir}\n\n")
            
            # 用 dictionary 手动统计每种 organism 出现次数。
            organisms = {}
            for info in seq_info:
                org = info['organism']
                organisms[org] = organisms.get(org, 0) + 1
            
            f.write("Sequences by Organism:\n")
            for org, count in sorted(organisms.items()):
                f.write(f"  {org}: {count}\n")
        
        # 重要：这里写出的 Total Sequences 是“通过筛选的 metadata 数量”，
        # 并不保证每一条 FASTA/GenBank 文件都实际下载成功。
        logger.info(f"Metadata saved to {self.output_dir}/")
    
    def create_combined_fasta(self) -> None:
        """依次读取 fasta/ 中的单个文件，拼成 all_rsv_genomes.fasta。"""
        logger.info("Creating combined FASTA file...")
        
        fasta_dir = os.path.join(self.output_dir, "fasta")
        combined_path = os.path.join(self.output_dir, "all_rsv_genomes.fasta")
        
        # 'w' 会重新建立 combined FASTA；再次运行会覆盖旧的 combined 文件。
        with open(combined_path, 'w') as combined_file:
            # os.listdir 的文件顺序不保证固定，因此合并后序列顺序也不保证固定。
            for filename in os.listdir(fasta_dir):
                if filename.endswith('.fasta'):
                    filepath = os.path.join(fasta_dir, filename)
                    with open(filepath, 'r') as individual_file:
                        combined_file.write(individual_file.read())
                        # 阅读提醒：上面 read() 后文件指针已经在末尾，下面第二次
                        # read() 会得到空字符串。因此这个条件实际上总为 True，
                        # 程序会在每个文件后加换行。通常不会破坏 FASTA，但写法有缺陷。
                        if not individual_file.read().endswith('\n'):
                            combined_file.write('\n')
        
        logger.info(f"Combined FASTA saved to {combined_path}")
    
    def analyze_dataset(self, seq_info: List[Dict]) -> Dict:
        """计算描述性统计；这不是序列比对、变异分析或系统发育分析。"""
        logger.info("Analyzing dataset...")
        
        if not seq_info:
            return {}
        
        # 建立最终会写入 dataset_analysis.json 的嵌套 dictionary。
        analysis = {
            'total_sequences': len(seq_info),
            'length_stats': {
                'min': min(int(info['length']) for info in seq_info),
                'max': max(int(info['length']) for info in seq_info),
                'mean': sum(int(info['length']) for info in seq_info) / len(seq_info)
            },
            'organisms': {},
            'years': {},
            'length_distribution': {}
        }
        
        # 统计每个 organism 标签的记录数。
        for info in seq_info:
            org = info['organism']
            analysis['organisms'][org] = analysis['organisms'].get(org, 0) + 1
        
        # 按 GenBank create_date 的年份统计“记录建立年份”。
        # 这不是病毒样本的 collection year。
        for info in seq_info:
            date_str = info.get('create_date', '')
            if date_str:
                try:
                    year = date_str.split('/')[0]  # 假设日期格式为 YYYY/MM/DD。
                    analysis['years'][year] = analysis['years'].get(year, 0) + 1
                except:
                    # 日期无法解析时静默跳过；裸 except 也可能隐藏其他意外错误。
                    pass
        
        # 把长度向下归入 100 nt 区间，例如 15,227 -> 15,200 bin。
        for info in seq_info:
            length = int(info['length'])
            length_bin = (length // 100) * 100  # Group by 100bp bins
            analysis['length_distribution'][length_bin] = analysis['length_distribution'].get(length_bin, 0) + 1
        
        return analysis
    
    def run_complete_download(self) -> None:
        """总控制器（orchestrator）：按照固定顺序运行整个 pipeline。"""
        logger.info("Starting complete RSV genome download process...")
        
        try:
            # Step 1：只搜索 accession。
            ids = self.search_rsv_genomes()
            if not ids:
                logger.error("No sequences found. Exiting.")
                return
            
            # Step 2：取得摘要 metadata，尚未下载核苷酸序列。
            seq_info = self.get_sequence_info(ids)
            if not seq_info:
                logger.error("Could not retrieve sequence information. Exiting.")
                return
            
            # Step 3：使用标题、organism 和长度做启发式筛选。
            filtered_info = self.filter_rsv_genomes(seq_info)
            if not filtered_info:
                logger.error("No RSV genomes found after filtering. Exiting.")
                return
            
            # Step 4：逐条下载 FASTA 与 GenBank 文件（最耗时的阶段）。
            self.download_sequences(filtered_info)
            
            # Step 5：保存“筛选后的 metadata”，不等于成功下载文件清单。
            self.save_metadata(filtered_info)
            
            # Step 6：只合并当前 fasta/ 文件夹中实际存在的文件。
            self.create_combined_fasta()
            
            # Step 7：根据 metadata 做简单汇总统计。
            analysis = self.analyze_dataset(filtered_info)
            
            # 把统计 dictionary 保存成 JSON。
            analysis_path = os.path.join(self.output_dir, "dataset_analysis.json")
            with open(analysis_path, 'w') as f:
                json.dump(analysis, f, indent=2)
            
            logger.info("Complete download process finished successfully!")
            logger.info(f"Results saved in: {self.output_dir}/")
            
        except Exception as e:
            # 未在内部处理的严重错误会在这里记录，然后重新抛出，使退出码表示失败。
            logger.error(f"Error in complete download process: {e}")
            raise


def main():
    """命令行入口：读取参数、检查依赖，并启动下载器。"""
    import argparse
    
    # argparse 把终端命令中的 --email、--output 等选项转换成 args 对象。
    parser = argparse.ArgumentParser(description="Download RSV whole genome sequences from GenBank")
    parser.add_argument("--email", required=True, help="Your email address (required by NCBI)")
    parser.add_argument("--output", default="rsv_genomes", help="Output directory (default: rsv_genomes)")
    parser.add_argument("--max-results", type=int, default=10000, help="Maximum number of search results")
    # 重要限制：args.max_results 在下面没有传给 search_rsv_genomes；
    # 因此本版本中 --max-results 虽然能被解析，实际上不会改变搜索上限。
    
    args = parser.parse_args()
    
    # 至少需要 Biopython 或 requests 中的一个；两者都没有时无法访问 NCBI。
    missing_deps = []
    if not BIOPYTHON_AVAILABLE and not REQUESTS_AVAILABLE:
        missing_deps.append("biopython OR requests")
    
    if missing_deps:
        print(f"Missing required dependencies: {', '.join(missing_deps)}")
        print("Please install with: pip install biopython requests pandas")
        sys.exit(1)
    
    # 建立 downloader 对象，然后启动完整 pipeline。
    downloader = RSVGenBankDownloader(args.email, args.output)
    downloader.run_complete_download()


# 只有“直接运行这个 .py 文件”时才执行 main()。
# 如果它被另一个 Python 文件 import，则不会自动开始下载。
if __name__ == "__main__":
    main()
