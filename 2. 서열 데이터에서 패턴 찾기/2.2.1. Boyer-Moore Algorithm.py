# 보이어-무어 클래스 정의
class BoyerMoore:
    
    # alphabet, pattern 객체 정의
    def __init__(self, alphabet, pattern):
        self.alphabet = alphabet
        self.pattern = pattern
        self.preprocess()
    
    # preprocess 함수 정의
    def preprocess(self):
        self.process_bcr()
        self.process_gsr()
        
    # 잘못된 문자 규칙
    def process_bcr(self):
        self.occ = {}
        for symb in self.alphabet:
            self.occ[symb] = -1
        for j in range(len(self.pattern)):
            c = self.pattern[j]
            self.occ[c] = j
    
    # 좋은 접미사 규칙
    def process_gsr(self):
        self.f = [0] * (len(self.pattern) + 1)
        self.s = [0] * (len(self.pattern) + 1)
        i = len(self.pattern)
        j = len(self.pattern)+1
        self.f[i] = j
        while i > 0:
            while j <= len(self.pattern) and self.pattern[i-1] != self.pattern[j-1]:
                if self.s[j] == 0: 
                    self.s[j] = j-i;
                j = self.f[j]
            i -= 1
            j -= 1
            self.f[i] = j
        j = self.f[0]
        for i in range(len(self.pattern)):
            if self.s[i] == 0:
                self.s[i] = j
            if i == j:
                j = self.f[j]
                
    # 주어진 패턴에 대한 서열 데이터 검색
    def search_pattern(self,text):
        res = []
        i = 0
        while i <= len(text) - len(self.pattern):
            j = len(self.pattern) - 1
            while j >= 0 and self.pattern[j] == text[j+i]:
                j -= 1
            if j < 0:
                res.append(i)
                i += self.s[0]
            else:
                c = text[j+i]
                i += max(self.s[j+1], j-self.occ[c])
        return res

# test( ) 함수 정의
def test():
    bm = BoyerMoore("ACTG", "ACCA")
    print(bm.search_pattern("ATAGAACCAATGAACCATGATGAACCATGGATACCCAACCACC"))

test()
