if __name__ == '__main__':
    if __package__ is None:
        import sys
        from os import path
        sys.path.append( path.dirname( path.dirname( path.abspath(__file__))))
        from matrix_module.matrixOperations import MOperations
        from galois_field.gfp import GFP
        from galois_field.simple import Simple
    else:
        from ..matrix_module.matrixOperations import MOperations
        from ..galois_field.gfp import GFP
        from galois_field.simple import Simple

import unittest

class TestMatrix(unittest.TestCase):

    def assertEqualWithInaccuaracy(self, matrixA, matrixB):
        row = len(matrixA)
        column = len(matrixA[0])
        for i in range(row):
            for j in range(column):
                if (abs(matrixA[i][j] - matrixB[i][j]) > 1e-14):
                    return False
        return True

    def test_mult_correct(self):
        operations = MOperations(Simple())
        matrixA = [[3, 2, 1],
                   [1, -7, 8],
                   [9, 1, 1]]
        matrixB = [[2, 1, 1],
                   [3, 2, -2],
                   [-1, -2, -3]]
        expectedMatrix = [[11, 5, -4],
                           [-27, -29, -9],
                           [20, 9, 4]]
        self.assertEqual(operations.matrixMult(matrixA, matrixB), expectedMatrix)

    def test_mult_error(self):
        operations = MOperations(Simple())
        matrixA = [[3, 2, 1],
                   [1, -7, 8],
                   [9, 1, 1]]
        matrixB = [[2, 1, 1],
                   [3, 2, -2]]
        with self.assertRaises(ArithmeticError):
            operations.matrixMult(matrixA, matrixB)

    # def test_solve_vandermond_rank4(self):
    #     operations = MOperations(Simple())
    #     vectorC = [[2], [-1], [-3], [-2]]
    #     vectorA = [[35], [2], [-50], [-13]]
    #     expectedResult = [[7], [4], [1], [2]]
    #     self.assertEqual(operations.solveVondermond(vectorA, vectorC), expectedResult)

    def test_solve_vandermond_rank3(self):
        operations = MOperations(Simple())
        vectorC = [[-3], [-7], [2]]
        vectorA = [[-16], [-108], [9]]
        expectedResult = [[11], [3], [-2]]
        self.assertEqual(operations.solveVondermond(vectorA, vectorC), expectedResult)

    def test_transpose(self):
        operations = MOperations(Simple())
        matrixA = [[3, 2, 1],
                   [1, -7, 8],
                   [9, 1, 1]]
        expectedMatrix = [[3, 1, 9],
                          [2, -7, 1],
                          [1, 8, 1]]
        self.assertEqual(operations.transposeMatrix(matrixA), expectedMatrix)

    def test_inverse(self):
        operations = MOperations(Simple())
        matrixA = [[7.0, 1.0, 2.0],
                   [-5.0, 2.0, 12.0],
                   [2.0, 3.0, 1.0]]
        identityMatrix = [[1.0, 0.0, 0.0],
                          [0.0, 1.0, 0.0],
                          [0.0, 0.0, 1.0]]
        inversedMatrix =  operations.getMatrixInverse(matrixA)
        calculatedIdentityMatrix  = operations.matrixMult(matrixA, inversedMatrix)
        self.assertTrue(self.assertEqualWithInaccuaracy(identityMatrix, calculatedIdentityMatrix))

    def test_generate_vondermond(self):
        operations = MOperations(Simple())
        vectorC = [[2], [-1], [-3], [-2]]
        expectedMatrix = [[1.0, 2.0, 4.0, 8.0],
                          [1.0, -1.0, 1.0, -1.0],
                          [1.0, -3.0, 9.0, -27.0],
                          [1.0, -2.0, 4.0, -8.0]]
        generatedMatrix = operations.generateVondermondMatrix(vectorC)
        self.assertTrue(self.assertEqualWithInaccuaracy(expectedMatrix, generatedMatrix))

    def test_solve_straight_way(self):
        operations = MOperations(Simple())
        vectorA = [[35.0], [2.0], [-50.0], [-13.0]]
        vondermondMatrix = [[1.0, 2.0, 4.0, 8.0],
                            [1.0, -1.0, 1.0, -1.0],
                            [1.0, -3.0, 9.0, -27.0],
                            [1.0, -2.0, 4.0, -8.0]]
        expectedResult = [[7.0], [4.0], [1.0], [2.0]]
        inversedMatrix = operations.getMatrixInverse(vondermondMatrix)
        calculatedResult = operations.matrixMult(inversedMatrix, vectorA)
        self.assertTrue(self.assertEqualWithInaccuaracy(expectedResult, calculatedResult))

    def test_compare_straight_and_vondermond(self):
        operations = MOperations(Simple())
        vectorA = [[35.0], [2.0], [-50.0], [-13.0]]
        vectorC = [[2.0], [-1.0], [-3.0], [-2.0]]
        vondermondResult = operations.solveVondermond(vectorA, vectorC)
        vondermondMatrix = operations.generateVondermondMatrix(vectorC)
        inversedMatrix = operations.getMatrixInverse(vondermondMatrix)
        straightResult = operations.matrixMult(inversedMatrix, vectorA)
        self.assertTrue(self.assertEqualWithInaccuaracy(vondermondResult, straightResult))

    def test_galue_field_prime_add(self):
        gfp = GFP(7)
        result = gfp.add(6, 5)
        self.assertEqual(result, 4)

    def test_galue_field_prime_sub(self):
        gfp = GFP(7)
        result = gfp.sub(5, 6)
        self.assertEqual(result, 6)

        gfp = GFP(2)
        result =  gfp.sub(0, 1)
        self.assertEqual(result, 1)

    def test_galue_field_prime_mult(self):
        gfp = GFP(7)
        result = gfp.mult(3, 4)
        self.assertEqual(result, 5)

        gfp = GFP(5)
        result = gfp.mult(4, 4)
        self.assertEqual(result, 1)

    def test_galue_field_prime_div(self):
        gfp = GFP(7)
        result = gfp.div(5, 6)
        self.assertEqual(result, 2)

        result = gfp.div(5, 1)
        self.assertEqual(result, 5)

    def test_galue_field_inverse(self):
        gfp = GFP(57)
        result = gfp.inverse(33)
        self.assertEqual(result, 7)

        gfp = GFP(5)
        result = gfp.inverse(3)
        self.assertEqual(result, 2)

    def test_galue_field_matrix_mult(self):
        gfp = GFP(7)
        operations = MOperations(gfp)
        matrixA = [[3, 2, 1],
                   [1, -7, 8],
                   [9, 1, 1]]
        matrixB = [[2, 1, 1],
                   [3, 2, -2],
                   [-1, -2, -3]]
        expectedMatrix = [[11, 5, -4],
                           [-27, -29, -9],
                           [20, 9, 4]]

        galueMatrixA = gfp.toGFP(matrixA)
        galueMatrixB = gfp.toGFP(matrixB)
        galueExpectedMatrix = gfp.toGFP(expectedMatrix)
        calcMatrix = operations.matrixMult(galueMatrixA, galueMatrixB)
        self.assertEqual(galueExpectedMatrix, calcMatrix)

    def test_galue_field_matrix_add(self):
        gfp = GFP(5)
        operations = MOperations(gfp)
        matrixA = [[6, 0, 5],
                   [12, -7, 8],
                   [3, 1, 1]]
        matrixB = [[2, 1, 1],
                   [3, 2, -2],
                   [-1, 0, 17]]
        galueExpectedMatrix = [[3, 1, 1],
                               [0, 0, 1],
                               [2, 1, 3]]
        galueMatrixA = gfp.toGFP(matrixA)
        galueMatrixB = gfp.toGFP(matrixB)
        calcMatrix = operations.matrixAddition(galueMatrixA, galueMatrixB)
        self.assertEqual(galueExpectedMatrix, calcMatrix)

    def test_galue_field_vondermond_rank4(self):
        gfp = GFP(7)
        operations = MOperations(gfp)
        vectorC = [[2], [-1], [-3], [-2]]
        vectorA = [[35], [2], [-50], [-13]]
        expectedResult = [[7], [4], [1], [2]]

        gVectorC = gfp.toGFP(vectorC)
        gVectorA = gfp.toGFP(vectorA)
        gExpResult = gfp.toGFP(expectedResult)
        self.assertEqual(operations.solveVondermond(gVectorA, gVectorC), gExpResult)

    def test_galue_field_vondermond_rank3(self):
        gfp = GFP(7)
        operations = MOperations(gfp)
        vectorC = [[-3], [-7], [2]]
        vectorA = [[-16], [-108], [9]]
        expectedResult = [[11], [3], [-2]]

        gVectorC = gfp.toGFP(vectorC)
        gVectorA = gfp.toGFP(vectorA)
        gExpResult = gfp.toGFP(expectedResult)
        self.assertEqual(operations.solveVondermond(gVectorA, gVectorC), gExpResult)

    def test_unique_vector_generator(self):
        operations = MOperations(Simple())
        vector = operations.createUniqueFloatVector(100, 3)

        vectorInt = operations.createUniqueIntegerVector(30, 20)
        print vectorInt

if __name__ == '__main__':
    unittest.main()
