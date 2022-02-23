class Max_finder:
    
    def max_beta(Br):
        maxBta = Br[0]
        i=0
        while i<len(Br):
            if len(maxBta[0]) < len(Br[i][0]):
                maxBta = Br[i]
            i=i+1
        return maxBta
    
    