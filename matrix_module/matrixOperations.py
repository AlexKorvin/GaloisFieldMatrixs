from copy import deepcopy
from random import *

class MOperations(object):
    ar = None

    def __init__(self, arithmetic):
        self.ar = arithmetic

    def readMatrixFromFile(self, fileName):
        matrix = []
        with open(fileName, 'r') as openFile:
            for line in openFile:
                line = line.strip() #remove extra spaces
                matrix.append(map(float, line.split())) #split off values and return them as list
        return matrix

    def readVectorFromFile(self, fileName):
        vector = []
        with open(fileName, 'r') as openFile:
            for line in openFile:
                line = line.strip()
                vector.extend(map(int, line.split()))
        return vector

    def printMatrixToFile(self, fileName, matrix):
        matrixToString = '\n'.join(' '.join(str(value) for value in raw) for raw in matrix)
        with open(fileName, 'w') as openFile:
            openFile.write(matrixToString)

    def printVectorToFile(self, fileName, vector):
        with open(fileName, 'w') as openFile:
            openFile.write(str(vector))

    def createUniqueFloatVector(self, size, valueBorder):
        vector = []
        tempSet = set()
        while (len(tempSet) < size):
            value = uniform(-valueBorder, valueBorder)
            if ((value not in tempSet) and (value != 0)):
                tempSet.add(value)
                vector.append([value])

        return vector

    def createUniqueIntegerVector(self, size, valueBorder):
        vector = []
        tempSet = set()
        while (len(tempSet) < size):
            value = randint(-valueBorder, valueBorder)
            if ((value not in tempSet) and (value != 0)):
                tempSet.add(value)
                vector.append([float(value)])

        return vector

    def transposeMatrix(self, matrix):
        return map(list, zip(*matrix))

    def assertRank(self, matrix):
        return false if len(matrix[0]) < 2 else true

    def assertSquareness(self, matrix):
        return false if len(matrix) != len(matrix[0]) else true

    def generateVondermondMatrix(self, vector):
        rank = len(vector)
        value = 1
        tempList = [value]
        vondermond = []

        for i in range(0, rank):
            for j in range(1, rank):
                value = self.ar.mult(value, vector[i][0])
                tempList.append(value)
            vondermond.append(tempList)
            value = 1
            tempList = [value]
        return vondermond

    def solveVondermond(self, vectA, vectC):
        shiftA = deepcopy(vectA)
        shiftX = []
        tempX = deepcopy(vectA)

        for iteration in range(0, len(vectA) - 1):
            for equation in range(iteration + 1, len(vectA)):
                top = self.ar.sub(shiftA[equation][0], shiftA[iteration][0])
                bot = self.ar.sub(vectC[equation][0], vectC[iteration][0])
                tempX[equation][0] = self.ar.div(top, bot)
            shiftA = deepcopy(tempX)

        shiftX = deepcopy(shiftA)
        for iteration in range(len(vectA) - 2, -1, -1):
            for equation in range(len(vectA) - 2, iteration - 1, -1):
                tempX[equation][0] = self.ar.sub(shiftX[equation][0], self.ar.mult(vectC[iteration][0], shiftX[equation + 1][0]))
            shiftX = deepcopy(tempX)
        return shiftX

    def matrixMult(self, matrix1, matrix2):
        multSumm = 0
        tempList = []
        resMatrix = []

        if (len(matrix2) != len(matrix1[0])):
            raise ArithmeticError("Matrix cannot be multiplied")
        else:
            rowNumber1 = len(matrix1)
            columnNumber1 = len(matrix1[0])
            rowNumber2 = len(matrix2)
            columnNumber2 = len(matrix2[0])

            for i in range(0, rowNumber1):
                for j in range(0, columnNumber2):
                    for z in range(0, columnNumber1):
                        mult = self.ar.mult(matrix1[i][z], matrix2[z][j])
                        multSumm = self.ar.add(multSumm, mult)
                    tempList.append(multSumm)
                    multSumm = 0
                resMatrix.append(tempList)
                tempList = []
            return resMatrix

    def matrixAddition(self, matrix1, matrix2):
        result = deepcopy(matrix1)

        if (len(matrix2) != len(matrix1)) or (len(matrix2[0]) != len(matrix1[0])):
            raise ArithmeticError("Matrix cannot be added")
        else:
            rowNumber1 = len(matrix1)
            columnNumber1 = len(matrix1[0])

            # rows
            for i in range(0, rowNumber1):
                # columns
                for j in range(0, columnNumber1):
                    result[i][j] = self.ar.add(matrix1[i][j], matrix2[i][j])
            return result

    def getMatrixMinor(self, matrix, i, j):
        return [row[:j] + row[j + 1:] for row in (matrix[:i] + matrix[i + 1:])]

    def getMatrixDeternminant(self, matrix):
        if len(matrix) == 2:
            #first = self.ar.mult(matrix[0][0], matrix[1][1])
            #second = self.ar.mult(matrix[0][1], matrix[1][0])
            #return self.ar.sub(first - second)
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

        determinant = 0
        for c in range(len(matrix)):
            minor = self.getMatrixMinor(matrix, 0, c)
            determinant += ((-1) ** c) * matrix[0][c] * self.getMatrixDeternminant(minor)
        return determinant

    def getMatrixInverse(self, matrix):
        determinant = self.getMatrixDeternminant(matrix)

        if len(matrix) == 2:
            return [[matrix[1][1] / determinant, -1 * matrix[0][1] / determinant],
                    [-1 * matrix[1][0] / determinant, matrix[0][0] / determinant]]

            #find matrix of cofactors
        cofactors = []
        for r in range(len(matrix)):
            cofactorRow = []
            for c in range(len(matrix)):
                minor = self.getMatrixMinor(matrix, r, c)
                cofactorRow.append(((-1) ** (r + c)) * self.getMatrixDeternminant(minor))
            cofactors.append(cofactorRow)
        cofactors = self.transposeMatrix(cofactors)
        for r in range(len(cofactors)):
            for c in range(len(cofactors)):
                cofactors[r][c] = cofactors[r][c] / determinant
        return cofactors
