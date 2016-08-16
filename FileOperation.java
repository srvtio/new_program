/*--- 最新版 2016.8.16 ---*/
/*
任意のファルダにこのプログラムを置き，実行すると，コード中で指定したフォルダに置いてある
フォルダの一覧を取得し，各フォルダのPorousDataRaw().datファイルから必要なデータを読み取って
１つのファイルにまとめるプログラム
/*
(--- 読み取りファイル ---)
各フォルダの
Parameter().dat
PorousDataRaw().dat

(--- 出力ファイル ---)
flux-porosity.dat ... 流量と空隙率の関係を出力するファイル
datalist.dat ... ほしいデータの一覧を出力するファイル

*/
import java.io.*;

public class FileOperation {

    public static void main(String[] args) {
	int iCellNumber1 = 32 ;
	int iCellNumber2 = 32 ;
	int seed         = 0;
	int itemp        = 0;
	int KnudsenPt    = 1 ;
	int iParallel    = 1 ;
	int iPattern     = 0;
	double l1 = 1.0;
	double rho1 = 0.00;
	double tau1 = 1.00;
	double grad_p = (rho1*tau1 - 1.0)/l1/2.0;
	
	double porosity      = 0.0;
	double Knudsen1      = 0.0;
	double Knudsen2      = 100.0;
	double dtemp         = 0.0;
	double FinalPorosity = 0.0;
	
	System.out.println(iPattern) ;

	double mom_p_tptal = 0.0;
	double sta_Fp = 0.0;
	double sta_Fm = 0.0;
	double por_area = 0.0;

	/*--- ファイルへの書き込み ---*/
	try{
	    File file1 = new File("datalist.dat") ;
	    PrintWriter pw1 = new PrintWriter(new BufferedWriter(new FileWriter(file1))) ;

	    File file2 = new File("flux-porosity.dat") ;
	    PrintWriter pw2 = new PrintWriter(new BufferedWriter(new FileWriter(file2))) ;

	    //ファイル名の一覧を取得する
	    File folder = new File("/home/kasahara/study/scatterer/random/8.10");
	    File files[] = folder.listFiles();

	    //取得した一覧を表示する
	    for (int i=0; i<files.length; i++) {
		System.out.println("ファイル" + (i+1) + "→" + files[i]);

		/*--- ファイルの読み取り ---*/
		try{
		    // ファイル0から読み取り
		    File readfile0 = new File(files[i]+"/Parameter(CN1_"+iCellNumber1+",CN2_"+iCellNumber2+").dat") ;
		    FileReader fr0 = new FileReader(readfile0) ;
		    StreamTokenizer st0 = new StreamTokenizer(fr0) ; 

		    st0.nextToken();
		    KnudsenPt = (int)st0.nval;
		    st0.nextToken();
		    iParallel = (int)st0.nval;
		    st0.nextToken();
		    iPattern = (int)st0.nval;
		    st0.nextToken();
		    porosity = st0.nval;

		    // Finalporosityの値を計算
		    seed = iPattern / 1000000;
		    itemp = iPattern - (seed*1000000 + 11600);
		    FinalPorosity = (double)itemp / 100.0;
		  
		    // ファイル1から読み取り
		    File readfile1 = new File(files[i]+"/PorousDataRaw(pt "+iPattern+",dpdx_"+grad_p+"0,poro_"+FinalPorosity+").dat") ;
		    FileReader fr1 = new FileReader(readfile1) ;
		    StreamTokenizer st1 = new StreamTokenizer(fr1) ; 

		    st1.nextToken() ;
		    Knudsen1 = st1.nval ;
		    st1.nextToken() ;
		    Knudsen2 = st1.nval ;
		    st1.nextToken() ;
		    double total_kp1_n = st1.nval ;
		    st1.nextToken() ;
		    double subporosity = st1.nval ;
		    st1.nextToken() ;
		    grad_p = st1.nval ;
		    st1.nextToken() ;
		    porosity = st1.nval ;
		    st1.nextToken() ;
		    double sta_F = st1.nval ;
		    st1.nextToken() ;
		    sta_Fp = st1.nval ;
		    st1.nextToken() ;
		    sta_Fm = st1.nval ;
		    st1.nextToken() ;
		    por_area = st1.nval ;

		    st1.nextToken() ;
		    double mom_l = st1.nval ;
		    st1.nextToken() ;
		    double mom_r = st1.nval ;
		    st1.nextToken() ;
		    double mom_t = st1.nval ;
		    st1.nextToken() ;
		    double mom_b = st1.nval ;
		    st1.nextToken() ;
		    double mom_total = st1.nval ;
		    st1.nextToken() ;
		    double p_l = st1.nval ;
		    st1.nextToken() ;
		    double p_r = st1.nval ;
		    st1.nextToken() ;
		    double p_sub = st1.nval ;
		    st1.nextToken() ;
		    double p_t = st1.nval ; 
		    st1.nextToken() ;
		    double p_b = st1.nval ;
		    st1.nextToken() ;
		    double p_total = st1.nval ;
		    st1.nextToken() ;
		    mom_p_tptal = st1.nval ;

		    st1.nextToken() ;
		    double rho_r = st1.nval ;
		    st1.nextToken() ;
		    double t_r = st1.nval ;

		    st1.nextToken() ;
		    double average_p_left = st1.nval ;
		    st1.nextToken() ;
		    double average_p_right = st1.nval ;
		    st1.nextToken() ;
		    double average_p_sub = st1.nval ;

		    fr1.close() ;

		}catch(Exception e){
		    System.out.println(e) ;
		}//読み取りの例外

		//データリストの作成
		pw1.printf("%8d  ", iPattern);
		pw1.printf("%8f  ", porosity);
		pw1.printf("%8f  ", mom_p_tptal); 
		pw1.printf("%8f  ", sta_Fp);
		pw1.printf("%8f  ", sta_Fm);
		pw1.printf("%8f %n", por_area);

		//空隙率と流量の関係		
		pw2.printf("%8f  ", porosity);
		pw2.printf("%8f  %n", -(sta_Fp+sta_Fm)/mom_p_tptal/2.0);

	    }//フォルダごとの処理

	    pw1.close();
	    pw2.close() ;

	}catch(IOException e){
	    System.out.println(e) ;
	}//書き込みの例外
    }
}
