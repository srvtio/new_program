/*--- 最新版 2016.8.1 ---*/
/*
2種類のKnudsen().datファイルから多孔質物体の配置の違いを検出して出力するプログラム
/*
(--- 読み取りファイル ---)
2種類のKnudsen().datファイル

(--- 出力ファイル ---)
compareConf.dat ... 異なる部分は1，同じ場合は値無し

*/
import java.io.*;

public class compareConf {

    public static void main(String[] args) {
	int incx,incy ;
	int iCellNumber1 = 32 ;
	int iCellNumber2 = 32 ;
	int iPattern1     = 1011670 ;
	int iPattern2     = 1011668 ;
	double x1,x2 ;
	double rho1 = 0.00 ;
	double tau1 = 1.00 ;
	double Knudsen1[][] = new double[iCellNumber1][iCellNumber2] ;
	double Knudsen2[][] = new double[iCellNumber1][iCellNumber2] ;

	/*--- ファイルの読み取り ---*/
	try{
	    File readfile1 = new File("KnudsenConf(Rho1_"+rho1+"0,Tau1_"+tau1+"0,pt "+iPattern1+").dat") ;
	    File readfile2 = new File("KnudsenConf(Rho1_"+rho1+"0,Tau1_"+tau1+"0,pt "+iPattern2+").dat") ;
	    FileReader fr1 = new FileReader(readfile1) ;
	    FileReader fr2 = new FileReader(readfile2) ;
	    StreamTokenizer st1 = new StreamTokenizer(fr1) ; 
	    StreamTokenizer st2 = new StreamTokenizer(fr2) ; 

	    File writefile1 = new File("compareConf.dat") ;
	    PrintWriter pw1 = new PrintWriter(new BufferedWriter(new FileWriter(writefile1))) ;

	     //多孔質配置のデータファイル読み取り
	    // nextToken()で次のデータを読み取る
	    for(incx=0 ; incx<iCellNumber1 ; incx++){
		for(incy=0 ; incy<iCellNumber2 ; incy++){
		    // ファイル1から読み取り
		    st1.nextToken() ;
		    x1 = st1.nval ;
		    st1.nextToken() ;
		    x2 = st1.nval ;
		    st1.nextToken() ;
		    Knudsen1[incx][incy] = st1.nval ;

		    // ファイル2から読み取り
		    st2.nextToken() ;
		    x1 = st2.nval ;
		    st2.nextToken() ;
		    x2 = st2.nval ;
		    st2.nextToken() ;
		    Knudsen2[incx][incy] = st2.nval ;

		    pw1.printf("%16f", x1) ;
		    pw1.printf("%16f", x2) ;
		    if(Knudsen1[incx][incy]!=Knudsen2[incx][incy]){
			pw1.printf("%16f %n", 1.0) ;
		    }else{
			pw1.printf("%n") ;
		    }
		}

		pw1.printf("%8s %n", "  ") ;
	    }

	    fr1.close() ;
	    fr2.close() ;

	    pw1.close() ;

	}catch(Exception e){
	    System.out.println(e) ;
	}
    }
}
