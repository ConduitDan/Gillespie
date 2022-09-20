from operator import eq
import random
import math
import re
from matplotlib import pyplot as plt

class StructuralFormula:
    def __init__(self,chemical,coefficient):
        self.chemical = chemical
        self.coefficient = coefficient
    def getChemical(self):
        return self.chemical
    def getCoefficient(self):
        return self.coefficient

class Reaction:
    def __init__(self,reactants = None,products = None, rate = None,equation = None):
        # Initalization, try to parse a full 'A + B -> C' equation
        if equation is not None:
            self.parseFullRxn(equation)
        else:
            # otherwise is reactans and products are string parse thoses individually
            if isinstance(reactants, str):
                self.reactants = self.parseRxn(reactants)
            else: 
                #other wise just assume we've been given structual formulas
                self.reactants = reactants

            if isinstance(products, str):
                self.products = self.parseRxn(products)
            else:
                self.products = products

        self.rate = rate
    def parseFullRxn(self,equation):
        (reactants,products) = equation.strip().split('->')
        self.reactants = self.parseRxn(reactants)
        self.products = self.parseRxn(products)
    def parseRxn(self,eqn):
        #first split by + signs
        chemicals = eqn.strip().split('+')
        structuralFormulas = []
        for chemical in chemicals:
            #find the (possible) coeffeient and the label
            # look for any number of digits grouped together
            coeffMatch = re.search('(\d+)',chemical)
            coeff = 1 # default to one
            # if we found a match
            if coeffMatch:
                # set the coefficent equal to it
                coeff = int(coeffMatch.group(1))
            label = ''

            labelMatch = re.search('([A-Za-z]+)',chemical)
            if labelMatch:
                label = labelMatch.group(1)
            else:
                raise Exception("Could not parse chemical equation %s"%chemical)
            
            structuralFormulas.append(StructuralFormula(label,coeff))
        return structuralFormulas



    def getPropensity(self,chemicalList):
        a = self.rate
        for chemical in self.reactants:
            for i in range(chemical.getCoefficient()):
                a *= chemicalList[chemical.getChemical()]
        return a

    def performReaction(self,chemicalList):
        # repeated code... an opportunity for refactoring
        for chemical in self.reactants:
            chemicalList[chemical.getChemical()]-= chemical.getCoefficient()
        for chemical in self.products:
            chemicalList[chemical.getChemical()]+= chemical.getCoefficient()
        
    
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
        if isinstance(self.chemicalList, list): 
            return tuple(self.chemicalList)
        if isinstance(self.chemicalList,dict):
            return tuple(self.chemicalList.values())

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
        plt.plot(self.T,self.chems)        
        plt.savefig(self.fileName+".png")
        plt.show()

if __name__=="__main__":
    # A + B -> C
    # rx0 = Reaction([StructuralFormula(0,1),StructuralFormula(1,1)],[StructuralFormula(2,1)],1)
    # C -> A + B
    # rx1 = Reaction([StructuralFormula(2,1)],[StructuralFormula(0,1),StructuralFormula(1,1)],1)
    # myGillespie = Gillespie([50,50,50],[rx0,rx1])

    rx0 = Reaction(equation='A + B -> C',rate = 1)
    rx1 = Reaction(equation='C->A+B',rate = 1)
    myGillespie = Gillespie({'A':50,'B':50,'C':50},[rx0,rx1])

    myObs = GillespieObserver("data")

    myGillespie.setObserver(myObs)

    myGillespie.run(10)
    myObs.Plot()
