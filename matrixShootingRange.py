from matrix_module import matrixOperations

def main():
    matrixA = matrixOperations.readMatrixFromFile('matrixA')
    matrixB = matrixOperations.readMatrixFromFile('file')
    #matrixAInv = matrixOperations.getMatrixInverse(matrixA)
    print matrixA
    print matrixB

    #matrixBlaBla = matrixOperations.matrixMult(matrixA, matrixAInv)
    #print matrixBlaBla
    #vectX = matrixOperations.defeatVolandemort(vectA, vectC)
    #matrixT = operations.transposeMatrix(matrixC)
    #matrixOperations.printVectorToFile('vectX', vectX)
    #operations.printMatrix(matrixC)

    additionMatr = matrixOperations.matrixAddition(matrixA, matrixB)
    print additionMatr

main()
