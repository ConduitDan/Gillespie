import random
import math
from matplotlib import pyplot as plt

class StructuralFormula:
    def __init__(self,chemicalNumber,coefficient):
        self.chemicalNumber = chemicalNumber
        self.coefficient = coefficient
    def getNumber(self):
        return self.chemicalNumber
    def getCoefficient(self):
        return self.coefficient

class Reaction:
    def __init__(self,reactants,products, rate):
        self.reactants = reactants
        self.products = products
        self.rate = rate
    def getPropensity(self,chemicalList):
        a = self.rate
        for chemical in self.reactants:
            for i in range(chemical.getCoefficient()):
                a *= chemicalList[chemical.getNumber()]
        return a

    def performReaction(self,chemicalList):
        # repeated code... an opportunity for refactoring
        for chemical in self.reactants:
            chemicalList[chemical.getNumber()]-= chemical.getCoefficient()
        for chemical in self.products:
            chemicalList[chemical.getNumber()]+= chemical.getCoefficient()

class GillespieObserver:
    def __init__(self,fileName):
        self.fileName = fileName
        self.file  = open(fileName+".txt",'w')
        self.T = []
        self.chems = []
    def __del__(self):
        self.file.close()
    def takeMeasurement(self,gillespie):
        self.file.write(str(gillespie.getT())+" ")
        self.T.append(gillespie.getT())
        self.file.write(str(gillespie.getChemicals())+"\n")
        self.chems.append(gillespie.getChemicals())
    def Plot(self):
        fig,ax = plt.subplots()
        for i in range(len(self.chems[0])):
            #plt.plot(self.T,list(frame[i] for frame in self.chems))
            plt.plot(self.T,self.chems)        
        plt.savefig(self.fileName+".png")
        plt.show()



        
    
class Gillespie:
    def __init__(self,initialValues,reactions):
        self.initialValues = initialValues
        self.chemicalList = initialValues
        self.reactions = reactions
        self.propensities = [0.0] * len(reactions) 
        self.totalP = 0
        self.T = 0
        self.observer = None
    def setObserver(self,obs):
        self.observer = obs 
    def reset(self):
        self.T = 0
        self.chemicalList = self.initialValues
    
    def getT(self):
        return self.T
    
    def getChemicals(self):
        return tuple(self.chemicalList)

    def run(self,maxT):
        while self.T<maxT:
            self.step()

    def step(self):
        self.calculatePropensities()
        self.updateTime()
        RxNum = self.chooseReaction()
        self.reactions[RxNum].performReaction(self.chemicalList)
        if self.observer is not None:
            self.observer.takeMeasurement(self)
        
    def calculatePropensities(self):
        self.totalP = 0
        for i in range(len(self.reactions)):
            self.propensities[i] = self.reactions[i].getPropensity(self.chemicalList)
            self.totalP += self.propensities[i]
    
    def updateTime(self):
        dt = 1/self.totalP*math.log(1/random.random())
        self.T += dt
    
    def chooseReaction(self):
        r = random.random()
        RxNum = 0
        tally = 0
        for a in self.propensities:
            tally += a
            if tally/self.totalP>r:
                break
            RxNum+=1
        return RxNum

class ReactionLibrary:
    def SimpleEquilibrium(N):
        # A + B -> C
        rx0 = Reaction([StructuralFormula(0,1),StructuralFormula(1,1)],[StructuralFormula(2,1)],1)
        # C -> A + B
        rx1 = Reaction([StructuralFormula(2,1)],[StructuralFormula(0,1),StructuralFormula(1,1)],1)
        myGillespie = Gillespie([N,N,0],[rx0,rx1])
        return myGillespie

    def Brusselator(N,b,c):
        rxList = [Reaction([],[StructuralFormula(0,1)],N)]
        rxList.append(Reaction([StructuralFormula(0,1)],[],1))
        rxList.append(Reaction([StructuralFormula(0,1)],[StructuralFormula(1,1)],b))
        rxList.append(Reaction([StructuralFormula(0,2),StructuralFormula(1,1)],[StructuralFormula(0,3)],c/(N*N)))
        return  Gillespie([N,N],rxList)

if __name__=="__main__":

    myObs = GillespieObserver("brussData")
    myGillespie = ReactionLibrary.Brusselator(1000,2.2,1)
    myGillespie.setObserver(myObs)
    myGillespie.run(20)
    myObs.Plot()