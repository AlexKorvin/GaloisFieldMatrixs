from copy import deepcopy

class GFP(object):
    p = None

    def __init__(self, p):
        self.p = p

    def mult(self, a, b):
        resMult = (a * b) % self.p
        return resMult

    def add(self, a, b):
        resAdd = (a + b) % self.p
        return resAdd

    def sub(self, a, b):
        resSub = (a - b) % self.p
        return resSub

    def inverse(self, b):
        a = self.p
        xPr = 1
        yCur = 1
        xCur = 0
        yPr = 0
        temp = 0

        while True:
            temp = a % b
            q = a // b
            if temp is 0:
                return yCur
            a = b
            b = temp

            tempX = xPr - q * xCur
            xPr = xCur
            xCur = tempX

            tempY = yPr - q * yCur
            yPr = yCur
            yCur = tempY

    def div(self, a, b):
        divGalois = (a * self.inverse(b)) % self.p
        return divGalois

    def toGFP(self, matrix):
        result = deepcopy(matrix)
        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                result[i][j] = matrix[i][j] % self.p
        return result

    def __repr__(self):
        return "GFP(" + str(self.p) + ")"    
