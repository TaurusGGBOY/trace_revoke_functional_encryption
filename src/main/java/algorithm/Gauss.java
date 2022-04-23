package algorithm;

import it.unisa.dia.gas.jpbc.Element;
import it.unisa.dia.gas.jpbc.Pairing;

import java.util.*;

public class Gauss {

    Element[][] matrix;
    Pairing pairing;

    public Gauss(TraceRevokeFE.IdTree tree, Pairing pairing) {
        if (tree.gamma.isEmpty()) return;
        matrix = new Element[tree.gamma.size()][tree.gamma.get(0).theta.size()];
        for (int i = 0; i < tree.gamma.size(); i++) {
            for (int j = 0; j < tree.gamma.get(0).theta.size(); j++) {
                matrix[i][j] = tree.gamma.get(i).theta.get(j);
            }
        }
        this.pairing = pairing;
    }

    List<Element> getOneSolutionVector() {
        // 1. 处理为上三角矩阵
        transUpperTriangleMatrix();
        // 2. 递归求解第一行
        Element[] vec = solve(0);
        return new ArrayList<>(Arrays.asList(vec));
    }

    Element getMatrixGcd(int i) {
        Element gcd = pairing.getZr().newOneElement();
        for (int j = 0; j < matrix[i].length; j++) {
            if (matrix[i][j].equals(pairing.getZr().newZeroElement())) continue;
            gcd = gcd.mul(matrix[i][j]);
        }
        return gcd;
    }

    Element getVectorGcd(Element[] vec, int startCol, int endCol) {
        Element gcd = pairing.getZr().newOneElement();
        for (int i = startCol; i < endCol; i++) {
            if (vec[i].equals(pairing.getZr().newZeroElement())) continue;
            gcd = gcd.mul(vec[i]);
        }
        return gcd;
    }

    private Element[] solve(int i) {
        int col = matrix[0].length;
        int startCol = isEmptyLine(i);
        int nextStartCol = isEmptyLine(i + 1);
        Element[] nextVec = new Element[col];
        if (nextStartCol != col) nextVec = solve(i + 1);
        else Arrays.fill(nextVec, pairing.getZr().newZeroElement());

        // get sum of nextVec[j] * matrix[i][j]
        Element lastSum = pairing.getZr().newZeroElement();
        for (int j = nextStartCol; j < col; j++) {
            if (nextVec[j].equals(pairing.getZr().newZeroElement())) continue;
            lastSum = lastSum.add(nextVec[j].mul(matrix[i][j]));
        }

        // get gcd of the whole row
        Element gcd = getVectorGcd(matrix[i], startCol, nextStartCol);
        if (!lastSum.equals(pairing.getZr().newZeroElement())) {
            gcd = gcd.mul(lastSum);
        }
        Element sum = pairing.getZr().newZeroElement();

        // init vec to store the result vector
        Element[] vec = new Element[col];
        Arrays.fill(vec, pairing.getZr().newZeroElement());
        for (int j = startCol + 1; j < nextStartCol; j++) {
            if (nextVec[j].equals(pairing.getZr().newZeroElement())) continue;
            sum = sum.add(gcd);
            vec[j] = gcd.div(matrix[i][j]);
        }
        Element times = pairing.getZr().newOneElement();
        if (!lastSum.equals(pairing.getZr().newZeroElement())) {
            sum = sum.add(gcd);
            times = gcd.div(lastSum).negate().getImmutable();
        }
        for (int j = nextStartCol; j < col; j++) {
            if (nextVec[j].equals(pairing.getZr().newZeroElement())) continue;
            vec[j] = nextVec[j].mul(times).getImmutable();
        }
        // last compute the first one by others sum
        vec[startCol] = sum.div(matrix[i][startCol]).negate().getImmutable();
        return vec;
    }

    void rowExchange(int i, int j) {
        int col = matrix[0].length;
        Element[] temp = new Element[col];
        System.arraycopy(matrix[i], 0, temp, 0, col);
        System.arraycopy(matrix[j], 0, matrix[i], 0, col);
        System.arraycopy(temp, 0, matrix[i], 0, col);
    }

    void colExchange(int i, int j) {
        int row = matrix.length;
        for (int k = 0; k < row; k++) {
            Element temp = matrix[k][i];
            matrix[k][i] = matrix[k][j];
            matrix[k][j] = temp;
        }
    }

    void bubbleRow(int startRow) {
        int row = matrix.length;
        int minRow = startRow;
        int minCol = isEmptyLine(startRow);
        for (int i = startRow + 1; i < row; i++) {
            int temp = isEmptyLine(i);
            if (temp < minCol) {
                minRow = i;
                minCol = temp;
            }
        }
        rowExchange(startRow, minRow);
    }

    // find the first col index not equal to 0
    int isEmptyLine(int i) {
        if (i >= matrix.length) return matrix[0].length;
        int j = 0;
        for (; j < matrix[0].length; j++) {
            if (!matrix[i][j].equals(pairing.getZr().newZeroElement())) {
                break;
            }
        }
        return j;
    }

    int isEmptyVec(Element[] vec) {
        int j = 0;
        for (; j < vec.length; j++) {
            if (!vec[j].equals(pairing.getZr().newZeroElement())) {
                break;
            }
        }
        return j;
    }


    Element colGcd(int startRow, int startCol) {
        Element gcd = pairing.getZr().newOneElement();
        for (int j = startRow; j < matrix.length; j++) {
            if (matrix[j][startCol].equals(pairing.getZr().newZeroElement())) continue;
            gcd = gcd.mul(matrix[j][startCol]);
        }
        return gcd;
    }

    void updateRowByGcd(int startRow, int startCol, Element gcd) {
        for (int i = startRow; i < matrix.length; i++) {
            if (matrix[i][startCol].equals(pairing.getZr().newZeroElement())) continue;
            Element times = gcd.div(matrix[i][startCol]).getImmutable();
            for (int j = startCol + 1; j < matrix[0].length; j++) {
                matrix[i][j] = matrix[i][j].mul(times).getImmutable();
            }
        }
    }

    void transUpperTriangleMatrix() {
        for (int i = 0; i < matrix.length; i++) {
            bubbleRow(i);
            int startCol = isEmptyLine(i);
            Element gcd = colGcd(i, startCol).getImmutable();
            updateRowByGcd(i, startCol, gcd);
        }
    }
}