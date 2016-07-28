/*--- 最新版 2016.7.28 ---*/
// import java.io.File ;
// import java.io.FileWriter ;
// import java.io.BufferedWriter ;
// import java.io.PrintWriter ;
// import java.io.IOException ;
import java.io.*;

public class RandomObj {

    public static void main(String[] args) {
	int incx,incy;
	int iCellNumber1 = 32 ;
	int iCellNumber2 = 32 ;
	int iTotalCell   = iCellNumber1 * iCellNumber2 ;
	int seed         = 1 ;
	int itemp        = 0 ;
	int Number       = 0 ;
	int KnudsenPt    = 1 ;
	int iParallel    = 7 ;
	int iPattern     = 1011660 ;
	double porosity      = 0.0 ;
	double Knudsen1      = 0.0 ;
	double Knudsen2      = 100.0 ;
	double dtemp         = 0.0 ;
	double FinalPorosity = 0.50 ;

	int NumberCell[][] = new int[iCellNumber1][iCellNumber2] ;
	double Knudsen[][] = new double[iCellNumber1][iCellNumber2] ;

	/*--- Knudsen[][]の初期値設定 ---*/
	for(incx=0 ; incx<iCellNumber1 ; incx++){
	    for(incy=0 ; incy<iCellNumber2 ; incy++){
		Knudsen[incx][incy] = Knudsen2 ;
	    }
	}

	/*--- 空隙部分をKnudsen1で塗りつぶしていく繰り返し計算 ---*/
	for(int i=0 ; i<iTotalCell ; i++){

	    /*--- 空隙部分のセルに通し番号をつける ---*/
	    Number = 0 ;

	    for(incx=0 ; incx<iCellNumber1 ; incx++){
		for(incy=0 ; incy<iCellNumber2 ; incy++){
		
		    if(Knudsen[incx][incy] == Knudsen2){
			Number = Number + 1 ;
			NumberCell[incx][incy] = Number ;
		    }
		
		}
	    }

	    /*--- 空隙部分の中からランダムにセルを選んで，その部分のKnusenをKnudsen1にする ---*/
	    for(incx=0 ; incx<iCellNumber1 ; incx++){
		for(incy=0 ; incy<iCellNumber2 ; incy++){
		
		    dtemp = Math.ceil(Number * random(seed)) ;
		
		    if(NumberCell[incx][incy] == dtemp){
			Knudsen[incx][incy] = Knudsen1 ;
		    }
		
		}
	    }

	    porosity = (double)Number / (double)(iCellNumber1*iCellNumber2) ; 

	    if(porosity < FinalPorosity) break ;

	}

	System.out.println(porosity) ;

	// /*--- ファイルの読み取り ---*/
	// double a[] = new double[1000] ;
	// int i = 0 ;
	// try{
	//     File file = new File("Parameter(CN1_"+iCellNumber1+",CN2_"+iCellNumber2+").dat") ;
	//     FileReader fr = new FileReader(file) ; 
	//     StreamTokenizer st = new StreamTokenizer(fr) ; 

	//     // nextToken()で次の行を読み取る
	//     st.nextToken() ;
	//     KnudsenPt = (int)st.nval ;
	//     st.nextToken() ;
	//     iParallel = (int)st.nval ;
	//     st.nextToken() ;
	//     iPattern = (int)st.nval ;
	//     st.nextToken() ;
	//     porosity = st.nval ;

	//     fr.close() ;

	// }catch(Exception e){
	//     System.out.println(e) ;
	// }

	// System.out.println(KnudsenPt) ; 
	// System.out.println(iParallel) ; 
	// System.out.println(iPattern) ; 
	// System.out.println(porosity) ; 

	/*--- ファイルへの書き込み ---*/
	try{
	    File file1 = new File("Parameter(CN1_"+iCellNumber1+",CN2_"+iCellNumber2+").dat") ;
	    File file2
		= new File("KnudsenData(CN1_"+iCellNumber1+",CN2_"+iCellNumber2+",pt"+iPattern+").dat") ;
	    File file3 = new File("porous_conf.gp") ;

	    PrintWriter pw1 = new PrintWriter(new BufferedWriter(new FileWriter(file1))) ;
	    PrintWriter pw2 = new PrintWriter(new BufferedWriter(new FileWriter(file2))) ;
	    PrintWriter pw3 = new PrintWriter(new BufferedWriter(new FileWriter(file3))) ;

	    //パラメータファイルの作成
	    pw1.printf("%2d %n", KnudsenPt) ;
	    pw1.printf("%2d %n", iParallel) ;
	    pw1.printf("%8d %n", iPattern) ;
	    pw1.printf("%8f %n", porosity) ;

	    //多孔質配置のデータファイル作成
	    for(incx=0 ; incx<iCellNumber1 ; incx++){
		for(incy=0 ; incy<iCellNumber2 ; incy++){
		    pw2.printf("%8d", incx) ;
		    pw2.printf("%8d", incy) ;
		    pw2.printf("%16f %n", Knudsen[incx][incy]) ;			
		}
		pw2.printf("%8s %n", "  ") ;
	    }

	    //gnuplot作成ファイルの作成
	    pw3.printf("set terminal eps %n") ;
	    pw3.printf("set pm3d map %n") ;
	    pw3.printf("set pm3d corners2color c2 %n") ;
	    pw3.printf("set size ratio 1 %n") ;
	    pw3.printf("set output '"+iPattern+".eps' %n") ;
	    pw3.printf("unset colorbox %n") ;
	    pw3.printf("set palette rgbformulae 22,13,-31 %n") ;
	    pw3.printf("set xrange [0:31] %n") ;
	    pw3.printf("set yrange [0:31] %n") ;
	    pw3.printf("set xlabel 'x1' %n") ;
	    pw3.printf("set ylabel 'x2' %n") ;
	    pw3.printf("splot 'KnudsenData(CN1_"+iCellNumber1+",CN2_"+iCellNumber2+",pt"+iPattern+").dat'u 1:2:3  palette notitle %n") ;
	    pw3.printf("unset output %n") ;
	    pw3.printf("reset %n") ;

	    pw1.close() ;
	    pw2.close() ;
	    pw3.close() ;

	}catch(IOException e){
	    System.out.println(e) ;
	}
    }

    // /*--- ファイルがすでに存在するか，普通のファイル，書き込み可能かを判断するメソッド ---*/
    // private static boolean checkBeforeWritefile(File file){
    // 	if (file.exists()){
    // 	    if (file.isFile() && file.canWrite()){
    // 		return true ;
    // 	    }
    // 	}

    // 	return false ;
    // }

    /*--- 乱数生成用メソッド ---*/
    static int ma[] = new int[56] ;
    static int inext  = 0 ;
    static int inextp = 0 ;
    static int iff = 0 ;
    public static double random(int idum) {
	int mbig = 1000000000 ;
	int mseed = 161803398 ;
	int mz = 0 ;
	double fac = 1.0e-9 ;
	int  mj = 0 ;
	int mk,i,ii,k ;
	double rf = 0.0 ;

	if(idum < 0 || iff == 0) {
	    iff = 1 ;
	    mj = mseed - Math.abs(idum) ;
	    mj = mj % mbig ;
	    ma[55] = mj ;
	    mk = 1 ;

	    for(i = 1 ; i < 55 ; i++){
		ii = (21*i) % 55 ;
		ma[ii] = mk ;
		mk = mj - mk ;
		if(mk<mz) mk = mk + mbig ;
		mj = ma[ii] ;
	    }

	    for(k = 1 ; k < 5 ; k++) {
		for(i = 1 ; i < 56 ; i++){
		    ma[i] = ma[i] - ma[1+((i+30)%55)] ;
		    if(ma[i]<mz) ma[i] = ma[i] + mbig ;
		}
	    }

	    inext = 0 ;
	    inextp = 31 ;
	}

	while((rf < 1.0e-8) || (rf > 0.99999999)){
	    inext = inext + 1 ;
	    if(inext==56) inext = 1 ;
	    inextp = inextp + 1 ;
	    if(inextp == 56) inextp = 1 ;
	    mj = ma[inext] - ma[inextp] ;
	    if(mj < mz) mj = mj + mbig ;
	    ma[inext] = mj ;
	    rf = mj * fac ;
	}

	return rf ;

    }

}
