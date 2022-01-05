# 코돈을 아미노산으로 번역하기 위한 표준 유전 코드 딕셔너리 

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
    
    
# 생물학적 서열의 클래스

class MySeq:
  
    # MySeq 클래스의 요소 : seq, seq_type - default : DNA
    def __init__(self, seq, seq_type = "DNA"):                  
        self.seq = seq.upper()
        self.seq_type = seq_type
    
    # seq의 길이를 반환
    def __len__(self):                                           
        return len(self.seq)
    
    # seq의 n번쨰 요소를 반환
    def __getitem__(self, n):                                   
        return self.seq[n]
    
    # seq을 슬라이싱
    def __getslice__(self, i, j):                               
        return self.seq[i:j]
    
    # seq 반환
    def __str__(self):                                          
        return self.seq
    
    # seq의 seq_type 반환
    def get_seq_biotype(self):                                      
        return self.seq_type
    
    # seq의 정보 반환 - seq, biotype
    def show_info_seq(self):                                     
        print("Sequence: " + self.seq + " biotype: " + self.seq_type)
        
    # 서열 종류에 따른 허용문자 
    def alphabet(self):
        if(self.seq_type == "DNA"): return "ACGT"
        elif (self.seq_type == "RNA"): return "ACGU"
        elif (self.seq_type == "PROTEIN"): return "ACDEFGHIKLMNPQRSTVWY"
        else: return None

    # 서열 검증 
    def validate(self):
        alp = self.alphabet()                                   # alphabet() 메서드를 받음
        res = True
        i = 0
        while i < len(self.seq) and res:
            if self.seq[i] not in alp: res = False              # 서열이 허용문자 내에 없다면 res = False
            else: i += 1                                        # 있으면 계속 진행
        return res
    
    # DNA서열을 RNA서열로 바꿔주는 전사함수 
    def transcription(self):
        if(self.seq_type == "DNA"):
            return MySeq(self.seq.replace("T","U"), "RNA")      # seq_type을 RNA로 replace
        else:
            return None
    
    # DNA서열의 역상보서열을 구하는 함수 
    def reverse_comp(self):
        if(self.seq_type != "DNA"): return None
        comp = ""
        for c in self.seq:
            if (c == "A"): comp = "T" + comp
            elif (c == "T"): comp = "A" + comp
            elif (c == "C"): comp = "G" + comp
            elif (c == "G"): comp = "C" + comp
        return MySeq(comp, "DNA") 
    
  # 단백질을 만드는 번역함수 
    def translate(self, iniPos = 0):
        if(self.seq_type != "DNA"): return None
        seq_aa = ""
        for pos in range(iniPos, len(self.seq)-2, 3):
            cod = self.seq[pos:pos+3]
            seq_aa += translate_codon(cod)                        # MySeq 클래스의 외부함수 translate_codon() - 클래스 외부함수 접근 역시 가능함
        return MySeq(seq_aa, "PROTEIN")                           # seq_type : PROTEIN
      
 
