#!/usr/bin/env python
# coding: utf-8

# In[1]:


#1.DNA 서열이 유효한지 체크 p.109
def validate_dna(dna_seq):
    seqm = dna_seq.upper()
    valid = seqm.count("A") + seqm.count("T") + seqm.count("G") + seqm.count("C")
    if valid == len(seqm): return True
    else: return False
    
#2.서열에서 각 심볼의 빈도 계산 p.110
def frequency(seq):
    dic = {}
    for s in seq.upper():
        if s in dic: dic[s] += 1
        else: dic[s] = 1
    return dic    

#1.코돈을 아미노산으로 번역하기 위한 표준 유전 코드 딕셔너리 p.113
def translate_codon(cod):
    
    tc = {"GCT":"A", "GCC":"A", "GCG":"A",
          "TGT":"C", "TGC":"C",
          "GAT":"D", "GAC":"D",
          "GAA":"E", "GAG":"E",
          "TTT":"F", "TTC":"F",
          "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
          "CAT":"H", "CAC":"H",
          "ATA":"I", "ATT":"I", "ATC":"I",
          "AAA":"K", "AAG":"K",
          "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
          "ATG":"M", 
          "AAT":"N", "AAC":"N",
          "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
          "CAA":"Q", "CAG":"Q",
          "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
          "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
          "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
          "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
          "TGG":"W",
          "TAT":"Y", "TAC":"Y",
          "TAA":"_", "TAG":"_", "TGA":"_"}
          
    if cod in tc: return tc[cod]
    else: return None

#1.입력한 DNA 서열을 전사한 RNA 서열을 만드는 함수 p.112
def transcription(dna_seq):
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    return dna_seq.upper().replace("T","U")

#2.DNA 서열의 역상보서열 p.112
def reverse_complement(dna_seq):
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    comp = ""
    for c in dna_seq.upper():
        if c == "A":
            comp = "T" + comp
        elif c == "T":
            comp = "A" + comp
        elif c == "G":
            comp = "C" + comp
        elif c == "C":
            comp = "G" + comp
    return comp

#4.DNA 서열에서 G/C 뉴클레오타이드의 퍼센트 반환 p.111
def gc_content(dna_seq):
    gc_count = 0
    for s in dna_seq:
        if s in "GCgc": gc_count += 1
    return gc_count / len(dna_seq)

#2.DNA 서열을 아미노산 서열로 번역 p.114
def translate_seq(dna_seq, ini_pos = 0):
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    seq_aa = ""
    
    for pos in range(ini_pos, len(seqm)-2, 3):
        cod = seqm[pos:pos+3]
        seq_aa += translate_codon(cod)
    return seq_aa

#1.역상보서열을 포함한 여섯개의 리딩 프레임에서 DNA서열을 계산 p.116
def reading_frames(dna_seq):
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    res = []
    #DNA서열
    res.append(translate_seq(dna_seq,0))
    res.append(translate_seq(dna_seq,1))
    res.append(translate_seq(dna_seq,2))
    #역상보서열
    rc = reverse_complement(dna_seq)
    res.append(translate_seq(rc,0))
    res.append(translate_seq(rc,1))
    res.append(translate_seq(rc,2))
    return res

#3.주어진 아미노산을 암호화하고 있는 각 코돈의 비율을 DNA서열로 표현 p.115
def codon_usage(dna_seq, aa):
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    dic = {}
    total = 0
    for i in range(0, len(seqm)-2, 3):
        cod = seqm[i:i+3]
        if translate_codon(cod) == aa:
            if cod in dic:
                dic[cod] += 1
            else:
                dic[cod] = 1
            total += 1
    if total > 0:
        for k in dic:
            dic[k] /= total
    return (dic, total)

#2.아미노산 서열에서 메싸이오닌 아미노산을 고려하여 가능한 단백질 리스트 생성 p.117
def all_proteins_rf(aa_seq):
    aa_seq = aa_seq.upper()
    current_prot = []
    proteins = []
    
    for aa in aa_seq:
        if aa == "_":
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                    current_prot = []
        else:
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    
    return proteins

#3.모든 오픈 리딩프레임에서 가능한 단백질 계산 p.117
def all_orfs(dna_seq):
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    rfs = reading_frames(dna_seq)
    res = []
    for rf in rfs:
        prots = all_proteins_rf(rf)
        for p in prots:
            res.append(p)
    return res

#4.모든 오픈 리딩프레임에서 가능한 단백질 계산 / 최소 크기로 걸러 정렬하는 리스트 반환 p.118
#정렬 삽입 함수
def all_orfs_ord(dna_seq, minsize = 0):
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    rfs = reading_frames(dna_seq)
    res = []
    for rf in rfs:
        prots = all_proteins_rf(rf)
        for p in prots:
            if len(p) > minsize: insert_prot_ord(p,res)
    return res

#정렬함수
def insert_prot_ord(prot, list_prots):
    i = 0                                                                                 # i = 0 초기화
    while i < len(list_prots) and len(prot) < len(list_prots[i]):                         # 기존에 있던 단백질 크기가 더 크면
        i += 1                                                                            # 인덱스 값 i 에 1씩 더함 - 내림차순 정렬
    list_prots.insert(i, prot)                                                            # 정렬 / 삽입
    
#여러 줄로 구성된 텍스트 파일에서 서열을 읽음. p.119
def read_seq_from_file(filename):
    fh = open(filename, "r")
    lines = fh.readlines()
    seq=""
    for l in lines:
        seq += l.replace("\n", "")
    fh.close()
    return seq

#서열을 파일에 입력. p.119
def write_seq_to_file(seq, filename):
    fh = open(filename, "w")
    fh.write(seq)
    fh.close()
    return None
