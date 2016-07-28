/*--- 最新版 2016.7.28 ---*/
/*
(--- 読み取りファイル ---)
無し

(--- 出力ファイル ---)
Parameter().dat ... メイン計算に必要なパラメータを記述したファイル
KnudsenData().dat ... 多孔質の配置を決めるファイル
porous_conf.gp ... gnuplotで多孔質の配置図を確認するためのファイル

*/
import java.io.* ;
import java.awt.* ;
import java.awt.event.* ;
import javax.swing.* ;

public class RandomObj extends JFrame{

    private JLabel lb, lb2 ;
    private JPanel[] pn = new JPanel[3] ;
    private JCheckBox ch1, ch2, tmp ;
    private JTextField tf ;

    public static void main(String[] args) {
	int incx,incy ;
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

	RandomObj sm = new RandomObj() ;

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

    public RandomObj()
    {
	super("多孔質物体配置プログラム") ;
	
	lb = new JLabel("初期ファイル作成") ;
	ch1 = new JCheckBox("通常") ;
	ch2 = new JCheckBox("衝突確率100%") ;
	tf = new JTextField() ;
	lb2 = new JLabel("並列数") ;

	for(int i=0 ; i<pn.length ; i++){
	    pn[i] = new JPanel() ;
	}

	setLayout(new GridLayout(3, 1)) ;

	pn[0].add(lb) ;
	pn[1].add(ch1) ;
	pn[1].add(ch2) ;
	tf.setPreferredSize(new Dimension(100, 20)) ;
	pn[2].add(lb2) ;
	pn[2].add(tf) ;

	add(pn[0]) ;
	add(pn[1]) ;
	add(pn[2]) ;
	
	ch1.addItemListener(new SampleItemListener()) ;
	ch2.addItemListener(new SampleItemListener()) ;
	tf.addActionListener(new SampleActionListener()) ;

	addWindowListener(new SampleWindowListener()) ;

	setSize(200,200) ;
	setVisible(true) ;
    }

    class SampleItemListener implements ItemListener
    {
	public void itemStateChanged(ItemEvent e)
	{
	    if(e.getStateChange() == ItemEvent.SELECTED){
		tmp = (JCheckBox) e.getSource() ;
		lb.setText(tmp.getText() + "を選びました．") ;
		try{
		    File file = new File("porous_conf.gp") ;
		    PrintWriter pw
			= new PrintWriter(new BufferedWriter(new FileWriter(file))) ;
		    pw.printf("おためし") ;
		    pw.close() ;
		}catch(IOException e2){
		    System.out.println(e2) ;
		}
	    }
	    else if(e.getStateChange() == ItemEvent.DESELECTED){
		tmp = (JCheckBox) e.getSource() ;
		lb.setText(tmp.getText() + "をやめました．") ;
	    }  
	}
    }

    class SampleActionListener implements ActionListener
    {
	public void actionPerformed(ActionEvent e)
	{
	    JTextField tmp2 = (JTextField) e.getSource() ;
	    lb2.setText(tmp2.getText() + "ですね") ;
	}
    }

    class SampleWindowListener extends WindowAdapter
    {
	public void windowClosing(WindowEvent e)
	{
	    System.exit(0) ;
	}
    }

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
