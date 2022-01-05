# 패턴이 처음 나타나는 위치를 반환하거나 발생하지 않으면 -1을 반환
def search_first_occ(seq, pattern):
    found = False
    i = 0
    while i < len(seq)-len(pattern) and not found:
        j = 0
        while j < len(pattern) and pattern[j]==seq[i+j]:
            j = j + 1
        if j == len(pattern):
            found = True
        else:
            i += 1
    if found: return i
    else: return -1
        
seqDNA = "ATAGAATAGATAATAGTC"
print(search_first_occ(seqDNA, "GAAT"))
print(search_first_occ(seqDNA, "TATA"))


# 서열에서 패턴이 일치하는지 검사하고 위치 리스트 반환
def search_all_occurrences(seq, pattern):
    res = []
    for i in range(len(seq)-len(pattern)+1):
        j = 0
        while j < len(pattern) and pattern[j]==seq[i+j]:
            j = j + 1
        if j == len(pattern):
            res.append(i)
    return res

seqDNA = "ATAGAATAGATAATAGTC"
print(search_all_occurrences(seqDNA, "AAT"))


# 사용자가 입력한 서열 데이터에 적용
def test_pat_search():
    seq = input("Input sequence: ")
    pat = input("Input pattern: ")
    print(pat, "occrs in the following positions: ")
    print(search_all_occurrences(seq, pat))
    
test_pat_search()
