　＃＃＃　形状及び結果の表示プログラム disp.c について　＃＃＃

＊コンパイル法
     cc disp.c -o disp -lX11 -lm
     gcc disp.c -o disp -lX11 -lm　など．

＊本文の記述との相違点
　　(1) Element(), Contour(), Vector() などの引数はほとんど外部変数と
　　　　しました．
　　(2) -vew の入力は，Ginit() ではなく Input() で行っています．
　　(3) -gf ファイルの指定のみを行えば，図1.3の標準formatのデ－タを表
　　　　示出来ます．ただし，隠面消去は行いません．隠面消去を行うのは
　　　　-d 3bem のオプションをつけたときだけです．
　　(4) 図1.3の標準formatだけでは処理できないため，-d 3fem オプション
　　　　のときの媒質境界面の表示は行いません．

＊フロッピィ－ディスク内のデ－タ表示のための入力例
　　CHAP1  disp -gf sphere.dat -d 3bem -vew 1. 1. 1.
           disp -gf sphere.dat -sf sphere.pot r -vew 1. 1. 1.
           disp -gf sphere.dat -vf sphere.edy i -vew 1. 1. 1.
           disp -vf sphere.edy r -vew 1. 1. 1.
    CHAP2  disp -gf magnet.dat
           disp -gf magnet.dat -sf magnet.pot
    CHAP5  disp -gf fort.59 -d 3fem -vew 1. 1. 1.
           disp -vf fort.31 r -vew 1. 1. 1.
    CHAP6  disp -gf f17 -d 3bem -vew 1. 1. 1.
    CHAP7  disp -gf sample.dat -d 3fem -vew 1. 1. 1.

                                                                 以上．
                                                                 
