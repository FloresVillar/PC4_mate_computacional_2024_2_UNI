import java.text.DecimalFormat;
import java.util.regex.*;



public class Simplex {
    private static DecimalFormat rd=new DecimalFormat("#.##");
    private static String []s={"s1","s2","s3","z"};
    private static String[]x={"x1","x2","x3","x4","s1","s2","s3","r"};
    // Método para realizar el algoritmo Simplex
    //--------------------------------------------------------------------------
    public static void simplexBase(double[][] tabla) {
        printTablaSimplexBase(tabla);
        int m = tabla.length; // Número de restricciones
        int n = tabla[0].length - 1; // Número de columnas n/2 = numero de variables

        // Paso 3: Ejecutar el algoritmo hasta llegar a la solución óptima
        while (true) {
            // Paso 3: Identificar la columna pivote (con el valor más negativo en la fila z)
            int colPivote = -1;
            double minValue = 0;
           
            // Buscar la columna pivote en la fila z
            for (int j = 0; j < n; j++) {
                if (tabla[m - 1][j] < minValue) {
                    minValue = tabla[m - 1][j];
                    colPivote = j;
                }
            }
            // Si no hay columna con valor negativo en la fila z, hemos llegado a la solución óptima
            if (colPivote == -1) {
                break;
            }
            // Paso 4: Identificar la fila pivote (buscar la menor razón)
            int filaPivote = -1;
            double minRatio = Double.MAX_VALUE;

            for (int i = 0; i < m - 1; i++) {
                double rhs = tabla[i][n];   //ultima columna
                double column = tabla[i][colPivote];
                if (column > 0) {
                    double ratio = rhs / column;
                    System.out.print("r/"+s[i]+" :\t"+rd.format(rhs)+"/"+rd.format(column)+" = "+rd.format(ratio)+"\n");
                    if (ratio < minRatio) {
                        minRatio = ratio;
                        filaPivote = i;
                    }
                }
            }
            // Paso 5: Intercambiar la variable básica(el titulo)
            // En este caso, la fila pivote se convierte en la nueva fila básica y la columna pivote como la nueva variable básica
            // Paso 6: Hacer el pivote
            double pivot = tabla[filaPivote][colPivote];
            System.out.println("haciendo 1 al pivote"+pivot+"imprimiendo la fila pivote");
            for (int j = 0; j <=n ; j++) {
                tabla[filaPivote][j] /= pivot; // Dividir la fila pivote por el valor del pivote
                System.out.print(rd.format(tabla[filaPivote][j])+" ");
            }
            System.out.println();
            printTablaSimplexBase(tabla);
            // Paso 7: Realizar operaciones en el resto de las filas
            for (int i = 0; i < m; i++) {
                if (i != filaPivote) {
                    double factor = tabla[i][colPivote];
                    System.out.println("fila "+s[i]+"-"+rd.format(factor)+"*fila pivote "+s[filaPivote]);
                    for (int j = 0; j <= n ; j++) {
                        tabla[i][j] -= factor * tabla[filaPivote][j];
                    }
                }
            }
            s[filaPivote] = x[colPivote];
            // Paso 7: Imprimir la tabla después de cada iteración (opcional)
            System.out.println("Tabla después de la iteración:");
            printTablaSimplexBase(tabla);
        }

        // Paso 8: Imprimir la solución óptima
        System.out.println("Solución óptima encontrada:");
        printTablaSimplexBase(tabla);
    }
    // Método para imprimir la tabla simplex
    //---------------------------------------------------------------------------
    public static void printTablaSimplexBase(double[][] tabla) {
        System.out.print("\t");
        for(int i=0;i<tabla[0].length;i++){
            System.out.print(x[i]+"\t");
        }
        System.out.println();
        for (int i = 0; i < tabla.length; i++) {
            System.out.print(s[i]+"\t");
            for (int j = 0; j < tabla[i].length; j++) {
                System.out.print(String.format("%.2f", tabla[i][j]) + "\t");
            }
            System.out.println();
        }
        System.out.println();
    }
    //-------------------------------------------------------------------------------

    public static void imprimir(double[][]m){
		System.out.println();
        for(int i=0;i<m.length;i++){
			for(int j=0;j<m[0].length;j++){
				System.out.printf("%1.2f\t",m[i][j]);
			}
			System.out.println();
		}
        System.out.println();
	}
    //-----------------------------------------------------------------------------
    public static void imprimirM(double[][]m){
		System.out.println();
        for(int i=0;i<m.length;i++){
			for(int j=0;j<m[0].length;j++){
				System.out.printf("%1.0f\t",m[i][j]);
			}
			System.out.println();
		}
        System.out.println();
	}
    //-----------------------------------------------------------------------------
    public static void simplex(String []res,int nRes,int nIncog) {
        // se recibe el numero de incognitas(columnas) y restricciones(filas), 1:simplex,2:simplex-dual,3-M
        //tabla[i->restricciones+1(z)][j->incognitas]= tablaArgumento(ocurrencia de coeficientes)
        //if(simplex) m columnas + a la derecha y ↓
        //concatenar  nueva matriz diagonal unitaria(restricciones) y al final tambien los r
        //2x1 + 1x2 + 1x3 ≤ 2
        //1x1 + 2x2 + 3x3 ≤ 5
        //2x1 + 2x2 + 1x3 ≤ 6
        //Z -3x1 -x2 -3x3 + 0
        double[][] matriz=new double[nRes+1][nIncog];
        double[][] terminoIn= new double[nRes+1][1];
        // Suponiendo que res[i] contiene la cadena de entrada (como la ecuación)
        for (int i = 0; i < nRes + 1; i++) {
            Matcher matcherCoef = Pattern.compile("([+-]?\\d+)(?=\\s*[a-zA-Z]\\d*)").matcher(res[i]); // Para coeficientes de las incógnitas
            int j = 0;
            // Captura los coeficientes de las incógnitas
            while (matcherCoef.find()) {
                if (j < nIncog) {
                    matriz[i][j] = Integer.parseInt(matcherCoef.group());
                    j++;
                }
            }   
            // Aquí nos aseguramos de que el puntero se ha movido más allá de los coeficientes
            // Ahora vamos a capturar el término independiente (que siempre debe estar al final)
            Matcher matcherTermino = Pattern.compile("(?<=\\s*(≤|<|=))\\s*([+-]?\\d+|[+-]0)").matcher(res[i]);
            if (matcherTermino.find()) {
                terminoIn[i][0] = Integer.parseInt(matcherTermino.group(2));
            }
        }
        //crear la matriz para las variables basicas
        imprimir(matriz);
        imprimir(terminoIn);
        double [][] basicas=new double[nRes+1][nRes];
        for(int i=0;i<nRes+1;i++){
            for(int j=0;j<nRes;j++){
                if(i==j){
                    basicas[i][j]=1;
                }
            }
        }
        imprimir(basicas);
        //ahora concatenar
        double[][] tabla=new double[nRes+1][nIncog+nRes+1];
        for(int i=0;i<nRes+1;i++){
            System.arraycopy(matriz[i], 0, tabla[i], 0, matriz[i].length);
        }
        for(int i=0;i<nRes+1;i++){
            System.arraycopy(basicas[i],0 ,tabla[i] , nIncog, basicas[i].length);
        }
        for(int i=0;i<nRes+1;i++){
            System.arraycopy(terminoIn[i],0 ,tabla[i] , nIncog+nRes, terminoIn[i].length);
        }
        imprimir(tabla);
        //crear las cabeceras S(vertical) y X(horizontal)
        String[]S=new String[nRes+1];
        for(int i=0;i<nRes;i++){
            S[i]="s"+(i+1);
        }
        S[nRes]="z";
        String[]X=new String[nIncog+nRes+1];
        for(int i=0;i<nIncog;i++){
            X[i]="x"+(i+1);
        }
        for(int i=nIncog;i<nIncog+nRes;i++){
            X[i]="s"+(i-nIncog+1);
        }
        X[nIncog+nRes]="r";
        printTablaSimplex(tabla,S,X);
        int m = tabla.length; // filas
        int n = tabla[0].length - 1; // Número de columnas n/2 = numero de variables

        // Paso 3: Ejecutar el algoritmo hasta llegar a la solución óptima
        while (true) {
            // Paso 3: Identificar la columna pivote (con el valor más negativo en la fila z)
            int colPivote = -1;
            double minValue = 0;
            // Buscar la columna pivote en la fila z
            for (int j = 0; j < n; j++) {
                if (tabla[m - 1][j] < minValue) {
                    minValue = tabla[m - 1][j];
                    colPivote = j;
                }
            }
            // Si no hay columna con valor negativo en la fila z, hemos llegado a la solución óptima
            if (colPivote == -1) {
                break;
            }
            // Paso 4: Identificar la fila pivote (buscar la menor razón)
            int filaPivote = -1;
            double minRatio = Double.MAX_VALUE;

            for (int i = 0; i < m - 1; i++) {
                double rhs = tabla[i][n];   //ultima columna
                double column = tabla[i][colPivote];
                if (column > 0) {
                    double ratio = rhs / column;
                    System.out.print("r/"+S[i]+" :\t"+rd.format(rhs)+"/"+rd.format(column)+" = "+rd.format(ratio)+"\n");
                    if (ratio < minRatio) {
                        minRatio = ratio;
                        filaPivote = i;
                    }
                }
            }
            // Paso 5: Intercambiar la variable básica(el titulo)
            // Paso 6: Hacer el pivote
            double pivot = tabla[filaPivote][colPivote];
            System.out.println("haciendo 1 al pivote"+pivot+"imprimiendo la fila pivote");
            for (int j = 0; j <=n ; j++) {
                tabla[filaPivote][j] /= pivot; // Dividir la fila pivote por el valor del pivote
                System.out.print(rd.format(tabla[filaPivote][j])+" ");
            }
            System.out.println();
            printTablaSimplex(tabla,S,X);
            // Paso 7: Realizar operaciones en el resto de las filas
            for (int i = 0; i < m; i++) {
                if (i != filaPivote) {
                    double factor = tabla[i][colPivote];
                    System.out.println("fila "+S[i]+"-"+rd.format(factor)+"*fila pivote "+S[filaPivote]);
                    for (int j = 0; j <= n ; j++) {
                        tabla[i][j] -= factor * tabla[filaPivote][j];
                    }
                }
            }
            S[filaPivote] = X[colPivote];
            // Paso 7: Imprimir la tabla después de cada iteración (opcional)
            System.out.println("Tabla después de la iteración:");
            printTablaSimplex(tabla,S,X);
        }

        // Paso 8: Imprimir la solución óptima
        System.out.println("Solución óptima encontrada:");
        printTablaSimplex(tabla,S,X);
    }
    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------
    public static void printTablaSimplex(double[][] tabla,String[]S,String[]X) {
        System.out.print("\t");
        for(int i=0;i<tabla[0].length;i++){
            System.out.print(X[i]+"\t");
        }
        System.out.println();
        for (int i = 0; i < tabla.length; i++) {
            System.out.print(S[i]+"\t");
            for (int j = 0; j < tabla[i].length; j++) {
                System.out.print(String.format("%.2f", tabla[i][j]) + "\t");
            }
            System.out.println();
        }
        System.out.println();
    }
    //-----------------------------------------------------------------------------------------------------
    public static void printTablaSimplexM(double[][] tabla,String[]S,String[]X) {
        System.out.print("\t");
        for(int i=0;i<tabla[0].length;i++){
            System.out.print(X[i]+"\t");
        }
        System.out.println();
        for (int i = 0; i < tabla.length; i++) {
            System.out.print(S[i]+"\t");
            for (int j = 0; j < tabla[i].length; j++) {
                System.out.print(String.format("%1.1f", tabla[i][j]) + "\t");
            }
            System.out.println();
        }
        System.out.println();
    }
    //-------------------------------------------------------------------------
    //------------------------- dual-------------------------------------------------
    public static void simplexDual(String []res,int nRes,int nIncog) {
        // se recibe el numero de incognitas(columnas) y restricciones(filas), 1:simplex,2:simplex-dual,3-M
        //tabla[i->restricciones+1(z)][j->incognitas]= tablaArgumento(ocurrencia de coeficientes)
        //if(simplex) m columnas + a la derecha y ↓
        //concatenar  nueva matriz diagonal unitaria(restricciones) y al final tambien los r
        //"1x1 + 4x4 >= 3.5",
        // "1x1 + 2x2 >=2.5",
        // "z -3x1 -8x2 + 0"
        double[][] matriz=new double[nRes+1][nIncog];
        double[][] terminoIn= new double[nRes+1][1];
        // Suponiendo que res[i] contiene la cadena de entrada (como la ecuación)
        for (int i = 0; i < nRes + 1; i++) {
            Matcher matcherCoef = Pattern.compile("([+-]?\\d+)(?=\\s*[a-zA-Z]\\d*)").matcher(res[i]); // Para coeficientes de las incógnitas
            int j = 0;
            // Captura los coeficientes de las incógnitas
            while (matcherCoef.find()) {
                if (j < nIncog) {
                    matriz[i][j] = Integer.parseInt(matcherCoef.group());
                    j++;
                }
            }   
            // Aquí nos aseguramos de que el puntero se ha movido más allá de los coeficientes
            // Ahora vamos a capturar el término independiente (que siempre debe estar al final)
            Matcher matcherTermino = Pattern.compile("(?<=\\s*(≤|<|=|>=|>|≥))\\s*([+-]?\\d+\\.\\d+|[+-]?\\d+|[+-]0)").matcher(res[i]);
            if (matcherTermino.find()) {
                terminoIn[i][0] = Double.parseDouble(matcherTermino.group(2));
            }
        }
        //crear la matriz para las variables basicas
        imprimir(matriz);
        imprimir(terminoIn);
        double[][]matrizDual=new double[nIncog][nRes];
        for(int i=0;i<matriz.length-1;i++){
            for(int j=0;j<matriz[0].length;j++){
                matrizDual[j][i] = matriz[i][j];
            }
        }
        imprimir(matrizDual);
        double [][] basicasDual=new double[nIncog][nIncog];
        for(int i=0;i<nIncog;i++){
            for(int j=0;j<nIncog;j++){
                if(i==j){
                    basicasDual[i][j]=1;
                }
            }
        }
        imprimir(basicasDual);
        //terminoInDual
        double[][]terminoInDual=new double[nIncog][1];
        for(int i=0;i<nIncog;i++){
            terminoInDual[i][0]=-matriz[nRes][i];
        }
        //ultima fila para la nueva tabla -terminoInd
        int l = matrizDual[0].length + basicasDual[0].length+1;
        double[][]ultimaFila = new double[1][l];
        for(int i=0;i<terminoIn.length;i++){
            ultimaFila[0][i]=-terminoIn[i][0];
        }
        imprimir(ultimaFila);
        //ahora concatenar
        nRes= matrizDual.length;
        nIncog = matrizDual[0].length;
        double[][] tabla=new double[nRes+1][nIncog+nRes+1];
        for(int i=0;i<nRes;i++){
            System.arraycopy(matrizDual[i], 0, tabla[i], 0, matrizDual[i].length);
        }
        for(int i=0;i<nRes;i++){
            System.arraycopy(basicasDual[i],0 ,tabla[i] , nIncog, basicasDual[i].length);
        }
        for(int i=0;i<nRes;i++){
            System.arraycopy(terminoInDual[i],0 ,tabla[i] , nIncog+nRes, terminoIn[i].length);
        }
        System.arraycopy(ultimaFila[0], 0, tabla[nRes],0,l);
        imprimir(tabla);
        //crear las cabeceras S(vertical) y X(horizontal)
        String[]S=new String[nRes+1];
        for(int i=0;i<nRes;i++){
            S[i]="s"+(i+1);
        }
        S[nRes]="z";
        String[]X=new String[nIncog+nRes+1];
        for(int i=0;i<nIncog;i++){
            X[i]="y"+(i+1);
        }
        for(int i=nIncog;i<nIncog+nRes;i++){
            X[i]="s"+(i-nIncog+1);
        }
        X[nIncog+nRes]="r";
        printTablaSimplex(tabla,S,X);
        int m = tabla.length; // filas
        int n = tabla[0].length - 1; // Número de columnas n/2 = numero de variables

        // Paso 3: Ejecutar el algoritmo hasta llegar a la solución óptima
        while (true) {
            // Paso 3: Identificar la columna pivote (con el valor más negativo en la fila z)
            int colPivote = -1;
            double minValue = 0;
           
            // Buscar la columna pivote en la fila z
            for (int j = 0; j < n; j++) {
                if (tabla[m - 1][j] < minValue) {
                    minValue = tabla[m - 1][j];
                    colPivote = j;
                }
            }
            // Si no hay columna con valor negativo en la fila z, hemos llegado a la solución óptima
            if (colPivote == -1) {
                break;
            }
            // Paso 4: Identificar la fila pivote (buscar la menor razón)
            int filaPivote = -1;
            double minRatio = Double.MAX_VALUE;

            for (int i = 0; i < m - 1; i++) {
                double rhs = tabla[i][n];   //ultima columna
                double column = tabla[i][colPivote];
                if (column > 0) {
                    double ratio = rhs / column;
                    System.out.print("r/"+S[i]+" :\t"+rd.format(rhs)+"/"+rd.format(column)+" = "+rd.format(ratio)+"\n");
                    if (ratio < minRatio) {
                        minRatio = ratio;
                        filaPivote = i;
                    }
                }
            }
            // Paso 5: Intercambiar la variable básica(el titulo)
            // En este caso, la fila pivote se convierte en la nueva fila básica y la columna pivote como la nueva variable básica
            // Paso 6: Hacer el pivote
            double pivot = tabla[filaPivote][colPivote];
            System.out.println("haciendo 1 al pivote"+pivot+"imprimiendo la fila pivote");
            for (int j = 0; j <=n ; j++) {
                tabla[filaPivote][j] /= pivot; // Dividir la fila pivote por el valor del pivote
                System.out.print(rd.format(tabla[filaPivote][j])+" ");
            }
            System.out.println();
            printTablaSimplex(tabla,S,X);
            // Paso 7: Realizar operaciones en el resto de las filas
            for (int i = 0; i < m; i++) {
                if (i != filaPivote) {
                    double factor = tabla[i][colPivote];
                    System.out.println("fila "+S[i]+"-"+rd.format(factor)+"*fila pivote "+S[filaPivote]);
                    for (int j = 0; j <= n ; j++) {
                        tabla[i][j] -= factor * tabla[filaPivote][j];
                    }
                }
            }
            S[filaPivote] = X[colPivote];
            // Paso 7: Imprimir la tabla después de cada iteración (opcional)
            System.out.println("Tabla después de la iteración:");
            printTablaSimplex(tabla,S,X);
        }

        // Paso 8: Imprimir la solución óptima
        System.out.println("Solución óptima encontrada:");
        printTablaSimplex(tabla,S,X);
    }
    //-------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    public static void simplexM(String []res,int nMenorIgual,int nMayorIgual,int nIncog,int nRes) {
        // se recibe el numero de incognitas(columnas) y restricciones(filas), 1:simplex,2:simplex-dual,3-M
        //tabla[i->restricciones+1(z)][j->incognitas]= tablaArgumento(ocurrencia de coeficientes)
        //if(simplex) m columnas + a la derecha y ↓
        //concatenar  nueva matriz diagonal unitaria(restricciones) y al final tambien los r
        //"2x1 + 1x2 ≤ 20",
        //"1x1 + 1x2 ≤ 18"
        //,"1x1 + 2x2 ≥ 12",
        //"z -5x1 -4x2 + 0"}
        double[][] matriz= new double[nRes+1][nIncog];
        double [][] basicas=new double[nRes+1][nRes];
        double [][] A=new double[nRes+1][nRes];
        double [][] terminoIn= new double[nRes+1][1];
        for(int i=0;i<nRes+1;i++){
            for(int j=0;j<nRes;j++){
                if(i==j){
                    basicas[i][j]=1;
                }
            }
        }
        for(int i=0;i<nRes+1;i++){
            for(int j=0;j<nRes;j++){
                if(i==j){
                    A[i][j]=1;
                }
            }
        }
        // Suponiendo que res[i] contiene la cadena de entrada (como la ecuación)
        for (int i = 0; i < nRes+1; i++) {
            Matcher matcherCoef = Pattern.compile("([+-]?\\d+)(?=\\s*[a-zA-Z]\\d*)").matcher(res[i]); // Para coeficientes de las incógnitas
            int j = 0;
            // Captura los coeficientes de las incógnitas
            while (matcherCoef.find()) {
                if (j < nIncog) {
                    matriz[i][j] = Integer.parseInt(matcherCoef.group());
                    j++;
                }
            }   
            // Aquí nos aseguramos de que el puntero se ha movido más allá de los coeficientes
            // Ahora vamos a capturar el término independiente (que siempre debe estar al final)
            Matcher matcherMenorIgual = Pattern.compile("(?<=\\s*(≤|<|=))\\s*([+-]?\\d+(\\.\\d+)?)").matcher(res[i]);
            if (matcherMenorIgual.find()) {
                terminoIn[i][0] = Double.parseDouble(matcherMenorIgual.group(2));
                basicas[i][i]=1;
                A[i][i] = 0;
                continue;
            }
            Matcher matcherMayorIgual = Pattern.compile("(?<=\\s*(>=|>|≥))\\s*([+-]?\\d+(\\.\\d+)?)").matcher(res[i]);
            if(matcherMayorIgual.find()){
                terminoIn[i][0] = Double.parseDouble(matcherMayorIgual.group(2));
                basicas[i][i]=-1;
                A[i][i]=1;
                A[nRes][i]=5000;//aqui va el M=5000
                continue;
            }
            if(res[i].matches(".*\\+\\s*0\\s*$")){
                terminoIn[i][0]=0;
            }
        }                                                       //"(?<=\\s*(≤|<|=|>=|>|≥))\\s*([+-]?\\d+\\.\\d+|[+-]?\\d+|[+-]0)"                                                                            
        //crear la matriz para las variables basicas
        imprimir(matriz);
        imprimir(basicas);
        imprimir(A);
        imprimir(terminoIn);
        //ahora concatenar
        double[][] tabla=new double[nRes+1][nIncog+2*nRes+1];
        for(int i=0;i<nRes+1;i++){
            System.arraycopy(matriz[i], 0, tabla[i], 0, matriz[i].length);
        }
        for(int i=0;i<nRes+1;i++){
            System.arraycopy(basicas[i],0 ,tabla[i] , nIncog, basicas[i].length);
        }
        for(int i=0;i<nRes+1;i++){
            System.arraycopy(A[i], 0, tabla[i], nIncog+nRes, A[i].length);
        }
        for(int i=0;i<nRes+1;i++){
            System.arraycopy(terminoIn[i],0 ,tabla[i] , nIncog+2*nRes, terminoIn[i].length);
        }
        imprimir(tabla);
        //crear las cabeceras S(vertical) y X(horizontal)
        String[]S=new String[nRes+1];
        for(int i=0;i<nRes;i++){
            if(A[i][i]>0){
                S[i]="A"+(i+1);
            }
            else{
                S[i]="s"+(i+1);
            }
        }
        S[nRes]="z";
        String[]X=new String[nIncog+2*nRes+1];
        for(int i=0;i<nIncog;i++){
            X[i]="x"+(i+1);
        }
        for(int i=nIncog;i<nIncog+nRes;i++){
            X[i]="s"+(i-nIncog+1);
        }
        for(int i=nIncog+nRes;i<nIncog+2*nRes;i++){
            X[i]="A"+(i-(nIncog+nRes)+1);
        }
        X[nIncog+2*nRes]="r";
        printTablaSimplexM(tabla,S,X);
        //hacer 0 el M (5000)
        for(int i=nIncog+nRes;i<nIncog+2*nRes;i++){
            if(tabla[nRes][i]>0){
                //operaciones elementales
                double factor= tabla[nRes][i];
                int ii=i -(nIncog+nRes);
                for(int j=0;j<tabla[0].length;j++){
                    tabla[nRes][j]=(tabla[nRes][j]-factor*tabla[ii][j]);
                }
            }
        }
        imprimirM(tabla);
        //el siguiente codigo proviene de simplexBase
        int m = tabla.length; // filas
        int n = tabla[0].length - 1; // Número de columnas n/2 = numero de variables

        // Paso 3: Ejecutar el algoritmo hasta llegar a la solución óptima
        while (true) {
            // Paso 3: Identificar la columna pivote (con el valor más negativo en la fila z)
            int colPivote = -1;
            double minValue = 0;
           
            // Buscar la columna pivote en la fila z
            for (int j = 0; j < n; j++) {
                if (tabla[m - 1][j] < minValue) {
                    minValue = tabla[m - 1][j];
                    colPivote = j;
                }
            }
            // Si no hay columna con valor negativo en la fila z, hemos llegado a la solución óptima
            if (colPivote == -1) {
                break;
            }
            // Paso 4: Identificar la fila pivote (buscar la menor razón)
            int filaPivote = -1;
            double minRatio = Double.MAX_VALUE;

            for (int i = 0; i < m - 1; i++) {
                double rhs = tabla[i][n];   //ultima columna
                double column = tabla[i][colPivote];
                if (column > 0) {
                    double ratio = rhs / column;
                    System.out.print("r/"+S[i]+" :\t"+rd.format(rhs)+"/"+rd.format(column)+" = "+rd.format(ratio)+"\n");
                    if (ratio < minRatio) {
                        minRatio = ratio;
                        filaPivote = i;
                    }
                }
            }
            // Paso 5: Intercambiar la variable básica(el titulo)
            // En este caso, la fila pivote se convierte en la nueva fila básica y la columna pivote como la nueva variable básica
            // Paso 6: Hacer el pivote
            double pivot = tabla[filaPivote][colPivote];
            System.out.println("haciendo 1 al pivote"+rd.format(pivot)+"imprimiendo la fila pivote");
            for (int j = 0; j <=n ; j++) {
                tabla[filaPivote][j] /= pivot; // Dividir la fila pivote por el valor del pivote
                System.out.print(rd.format(tabla[filaPivote][j])+"\t");
            }
            System.out.println("\n");
            printTablaSimplexM(tabla,S,X);
            // Paso 7: Realizar operaciones en el resto de las filas
            for (int i = 0; i < m; i++) {
                if (i != filaPivote) {
                    double factor = tabla[i][colPivote];
                    System.out.println("fila "+S[i]+"-"+rd.format(factor)+"*fila pivote "+S[filaPivote]);
                    for (int j = 0; j <= n ; j++) {
                        tabla[i][j] -= factor * tabla[filaPivote][j];
                    }
                }
            }
            S[filaPivote] = X[colPivote];
            // Paso 7: Imprimir la tabla después de cada iteración (opcional)
            System.out.println("Tabla después de la iteración:");
            printTablaSimplexM(tabla,S,X);
        }

        // Paso 8: Imprimir la solución óptima
        System.out.println("Solución óptima encontrada:");
        printTablaSimplexM(tabla,S,X);
    }
    //-----------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------
    // Método principal
    //-------------------------------------------------------------------------------
    public static void main(String[] args) {
        // Definir el tableau inicial basado en el problema de ejemplo
        // Formato: coeficientes de las variables y el RHS (última columna es RHS)
        double[][] coeficientes = {
            {2,     1,      3,      4,      1,      0,      0,      200},
            {3,     2,      4,      2,      0,      1,      0,      300},
            {1,     1,      2,      3,      0,      0,      1,      150},
            {-50, -40,      -70,    -60,    0,      0,      0,      0} // Fila z
        };
        // Ejecutar el algoritmo Simplex
        //simplexBase(coeficientes);
        String [] Tarea1Simplex = {"2x1 + 1x2 + 1x3 ≤ 2",
                                    "1x1 + 2x2 + 3x3 ≤ 5",
                                    "2x1 + 2x2 + 1x3 ≤ 6",            
                                    "Z-3x1 -1x2 -3x3 +0"};
        //simplex(Tarea1Simplex,3,3);
        String[] tarea2SimplexMuebles ={"2x1 + 1x2 + 3x3 + 4x4 <=200",
                                        "3x1 + 2x2 + 4x3 + 2x4 <=300", 
                                        "1x1 + 1x2 + 2x3 + 3x4 <=150" ,
                                    "z -50x1 -40x2 -70x3 -60x4 + 0"};
        //simplex2(tarea2SimplexMuebles, 03, 04);
        String[] tarea1SimplexDual = {"1x1 + 4x4 >= 3.5",
                                        "1x1 + 2x2 >=2.5",
                                        "z -3x1 -8x2 + 0"};
        //simplexDual(tarea1SimplexDual, 02, 2);
        String[] tarea2SimplexDual ={"2x1 + 1x2 + 1x3 ≥ 40",
                                    "1x1 + 1x2 + 1x3 ≥ 60", 
                             "z -70x1 -40x2 -90x3 + 0"};
        //simplexDual(tarea2SimplexDual, 02, 03);
        String[] tarea3SimplexDual={"1x1 + 2x2 + 3x3 ≥ 5",
                                    "2x1 + 2x2 + 1x3 ≥ 6",
                                    "z -3x1  -4x2 -5x3 + 0"};
        simplexDual(tarea3SimplexDual, 2, 3);
        String [] tareaSimplexM={"2x1 + 1x2 ≤ 20",
                                "1x1 + 1x2 ≤ 18"
                                ,"1x1 + 2x2 ≥ 12",
                                "z -5x1 -4x2 + 0"};
        simplexM(tareaSimplexM,2,1,2,3);
    }
}
