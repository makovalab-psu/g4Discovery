## copied from https://github.com/AnimaTardeb/G4Hunter/blob/master/G4Hunter.py

def BaseScore(line):
        item, liste=0, []
        #calcule le item de chaque base et la stock dans liste
        while ( item < len(line)):
            #a la fin d une sequence il est possible d avoir des GGG dans se cas
            # on verifie si la secore+1<len(line) car il ya un deuxieme G 
            #et 
            if (item < len(line) and (line[item]=="G" or line[item]=="g")):
                liste.append(1)
                #print liste
                if(item+1< len(line) and (line[item+1]=="G" or line[item+1]=="g")):
                    liste[item]=2
                    liste.append(2)
                    if (item+2< len(line) and (line[item+2]=="G" or line[item+2]=="g")):
                        liste[item+1]=3
                        liste[item]=3
                        liste.append(3)
                        if (item+3< len(line) and (line[item+3]=="G" or line[item+3]=="g")):
                            liste[item]=4
                            liste[item+1]=4
                            liste[item+2]=4
                            liste.append(4)
                            item=item+1
                        item=item+1
                    item=item+1
                item=item+1
                while(item < len(line) and (line[item]=="G" or line[item]=="g")):
                        liste.append(4)
                        item=item+1
            
            elif (item < len(line) and line[item]!="G" and line[item]!="g" and line[item]!= "C" and line[item]!="c" ):
                        liste.append(0)
                        item=item+1
                
            elif(item < len(line) and (line[item]=="C" or line[item]=="c")):
                liste.append(-1)
                if(item+1< len(line) and (line[item+1]=="C" or line[item+1]=="c" )):
                    liste[item]=-2
                    liste.append(-2)
                    if (item+2< len(line) and (line[item+2]=="C" or line[item+2]=="c" )):
                        liste[item+1]=-3
                        liste[item]=-3
                        liste.append(-3)
                        if (item+3< len(line) and (line[item+3]=="C" or line[item+3]=="c"  )):
                            liste[item]=-4
                            liste[item+1]=-4
                            liste[item+2]=-4
                            liste.append(-4)
                            item=item+1
                        item=item+1   
                    item=item+1
                item=item+1
                while(item < len(line) and (line[item]=="C" or line[item]=="c")):
                    liste.append(-4)
                    item=item+1
            
            else:
                    item=item+1 #la fin du la ligne ou y a des entrers
        return line, liste

def CalScore(liste, k):
        Score_Liste=[]
        #calcule de la moynne des scores pour toutes les sequences - les derniers k bases
        for i in range (len (liste)-(k-1)):
            #print (len(liste)), i, k
            j,Sum=0,0
            while (j<k):
                #print j, i
                Sum=Sum+liste[i]
                j=j+1
                i=i+1
            Mean=Sum/float(k)
            Score_Liste.append(Mean) 
        return Score_Liste
