// 最新版 2016.7.27
import java.io.File;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.IOException;

public class RandomObj {

    public static void main(String[] args) {
	int iCellNumber1 = 32 ;
	int iCellNumber2 = 32 ;

	double Knudsen1 = 0.0 ;
	double Knudsen2 = 100.0 ;

	int seed = 36 ;

	int incx,incy ;
	int Number = 0 ;

	int NumberCell[][] = new int[iCellNumber1][iCellNumber2] ;
	double Knudsen[][] = new double[iCellNumber1][iCellNumber2] ;

	int itemp = 0 ;
	double dtemp = 0.0 ;

	double porosity = 0.0 ;
	double FinalPorosity = 0.62 ;

	for(incx=0 ; incx<iCellNumber1 ; incx++){
	    for(incy=0 ; incy<iCellNumber2 ; incy++){
		Knudsen[incx][incy] = Knudsen2 ;
	    }
	}

	for(int i=0 ; i<1000 ; i++){
	    // 空隙部分のセルに通し番号をつける
	    Number = 0 ;

	    for(incx=0 ; incx<iCellNumber1 ; incx++){
		for(incy=0 ; incy<iCellNumber2 ; incy++){
		
		    if(Knudsen[incx][incy] == Knudsen2){
			Number = Number + 1 ;
			NumberCell[incx][incy] = Number ;
		    }
		
		}
	    }

	    // 空隙部分の中からランダムにセルを選んで，その部分のKnusenをKnudsen1にする
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
	
	// ファイルへの書き込み
	try{
	    File file = new File("KnudsenData(CN1_32,CN2_32,pt1544).dat");

	    if (checkBeforeWritefile(file)){
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));

		for(incx=0 ; incx<iCellNumber1 ; incx++){
		    for(incy=0 ; incy<iCellNumber2 ; incy++){
			pw.printf("%8d", incx) ;
			pw.printf("%8d", incy) ;
			pw.printf("%16f %n", Knudsen[incx][incy]) ;			
		    }
		    pw.printf("%8s %n", "  ") ;
		}

		pw.close();

	    }else{
		System.out.println("ファイルに書き込めません");
	    }
	}catch(IOException e){
	    System.out.println(e);
	}
    }

    // ファイルがすでに存在するか，普通のファイル，書き込み可能かを判断するメソッド
    private static boolean checkBeforeWritefile(File file){
	if (file.exists()){
	    if (file.isFile() && file.canWrite()){
		return true;
	    }
	}

	return false;
    }

    // 乱数生成用メソッド
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
