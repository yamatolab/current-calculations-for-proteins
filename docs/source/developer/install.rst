================================
DEPRECATED, applies to CURP v1.1
================================

このページでは、開発者のためのインストールの手順を紹介する.
以下ではLinuxとMacでのインストール手順を述べる.

安定板だけを試したい方は `こちら <./install>`__.

Install on Mac
===============

Macはバージョン毎にアプリケーションの揃え方が異なるので注意すること。
ちなみに、各バージョンの名前は次のように使う。

*  Mac OS X 10.6.x (Snow Leopard) ⇒ SLeopard
*  Mac OS X 10.7.x (Lion) ⇒ Lion
*  OS X 10.8.x (Mountain Lion) ⇒ MLion

また、MacでCurpを動かすために必要となるアプリケーションは次の通り：

*  Xcode Command Line Tools >= 4.5.4
*  Python ≧2.7.3
*  OpenMPI ≧1.4.3
*  gfortran ≧
*  Mercurial ≧2.2.3

出来るだけ最新のバージョンを揃えておいた方が良い.

Install Xcode
---------------

Mac上での開発環境(gcc等)を整えるためにXcodeをインストールする必要がある。
SLeopardではXcodeをインストールすれば、自動的にCommand Line Toolsも
インストールされる。LionとMLionではComman Line ToolsはXcodeの標準ライブラリとは別になっているので、別途インストールする必要がある。

Xcodeを立ち上げてから、メニューで、

::

   Xcode > Preference > Downloads > Components

で Command Line Toolsをインストール出来る。

Install Python
---------------

1. Download Python

   以下からPython-2.7.xの最新バージョンをダウンロードする:

   http://www.python.org/download/

   の「Python 2.7.5 Mac OS X 64-bit/32-bit x86-64/i386 Installer」
   を選ぶこと。

2. Install Python

   ダウンロードしたdmgファイルをダブルクリックして、
   他のアプリケーションと同じようにインストールすれば良い。

Install OpenMPI
----------------

SLeopardまではOpenMPIが標準でインストールされていた。
しかし、LionとMLionではOpenMPIはインストールされていないので、
自分でインストールしなければならない。
SLeopardでは以下の作業は必要無い。

1. 以下のサイトからOpenMPIをインストールする。

   http://www.open-mpi.org/software/ompi/v1.6/

   のopenmpi-1.6.tar.gz をダウンロードする。
   最新版があればそれを使うことにする。
   以降、バージョン番号は最新のものに読み替えること。

2. openmpiのソースコードをダウンロードしたディレクトリに移動して、
   次のようにコマンドを実行する::

      $ tar xzvf openmpi-1.6.tar.gz
      $ cd openmpi-1.6
      $ ./configure --prefix=<directory_path_to_want_to_install>
      $ make
      $ make install

3. 環境変数を設定する。

   bashrcに以下の設定を追加する::

      OPENMPI_PATH=<installed_directory_path>
      export PATH=$OPENMPI_PATH/bin:$PATH
      export LD_LIBRARY_PATH=$OPENMPI_PATH/lib:$LD_LIBRARY_PATH
      export DYLD_LIBRARY_PATH=$OPENMPI_PATH/lib:$DYLD_LIBRARY_PATH
      export MANPATH=$OPENMPI_PATH/share/man:$OPENMPI_PATH

Install GFortran
-----------------

1. gfortranを以下のURLからダウンロードする。
   64bit版を選ぶこと。

   http://gcc.gnu.org/wiki/GFortranBinaries

2. gfortranを解凍して、任意の場所に配置する。

   例：<DIRPATH_TO_GFORTRAN>

3. gfortranを使うための環境変数を設定する。

   /usr/local等に配置した場合には、この手順は必要無い。
   二つの環境変数を~/.bashrc等に追記する::

      GFORTRAN_PATH=<DIRPATH_TO_GFORTRAN>
      export PATH=$GFORTRAN_PATH/bin:$PATH
      export LD_LIBRARY_PATH=$GFORTRAN_PATH/lib:$LD_LIBRARY_PATH
      export DYLD_LIBRARY_PATH=$GFORTRAN_PATH/lib:$DYLD_LIBRARY_PATH

Install Mercurial
------------------

次のサイトからダウンロードしてインストールする．

http://mercurial.selenic.com/

Install CURP
-------------

1. ターミナル上で、Curpを置きたいディレクトリを作成して移動する。
   仮にインストールしたい場所を~/programs/curpにしたいのであれば、

   ::

      $ mkdir ~/programs/curp
      $ cd ~/programs/curp

2. curpのcloneをダウンロードする。

   1. で mkdir を行ったディレクトリで、次のようにコマンドを打つ。

   ::

      $ hg clone ssh://hg@bitbucket.org/takayamato/curp/wiki <1. で作ったディレクトリ>
   
..   次のようにコマンドを打つ。<user_name>は自分のBitBucketのアカウントである.
   ::

      $ hg clone https://<user_name>@bitbucket.org/takayamato/curp  ~/programs/curp
      警告: bitbucket.org の証明書 (fingerprint は 24:9c:45:8b:9c:aa:ba:55:4e:01
      :6d:58:ff:e4:28:7d:2a:14:ae:3b) 検証を省略(設定ファイルの hostfingerprints
      ないし web.cacerts 設定を確認してください)
      HTTP 認証を要求しました
         認証領域: Bitbucket.org HTTP
         ユーザ: <user_name>
         パスワード: 

   すると、上のように聞かれるので、パスワードをタイプすること。
   ダウンロードには時間が掛かる可能性があるので気長に待つ。

   もし, SSHの公開鍵を自分のアカウントに追加している場合にはパスワードを
   聞かれることはない．

3. CURP_HOMEを$HOME/.bashrc等に追記する。

   例：export CURP_HOME=$HOME/opt/curp

4. CURPを構築する。

   次のようにする。少々時間が掛かる。
   ::

      $ cd $CURP_HOME
      $ make

   CURPにはvirtualenvによるPython環境が必要であるが、
   これはPython環境の構築からCURPライブラリのビルドまでを全て行っている。

   ここで、 `make` に他のターゲットを与えることが出来る。

   fortranライブラリのコンパイルに使われるものはデフォルトではgfortranである。
   ``intel`` を指定するとifortが使われる::

      $ make intel

   また、並列計算バージョンを使いたくない時には、

   ::

      $make serial

   とすれば、pythonのmpiライブラリをビルドしない。

   環境を全て真っ新にしたい場合には、

   ::

      $ make clean

   とすれば良い。
      

6. (オプション)CURPのAPIドキュメントをビルドする。

   CURPの技術仕様ドキュメントを自動生成することが出来る::

      $ cd $CURP_HOME
      $ make apihtml

         # ここで、生成が完了する。

      $ cd $CURP_HOME/docs/_api/
      $ open index.html

   PDFファイルも生成することが出来る::

      $ cd $CURP_HOME
      $ make apipdf

         # ここで、生成が完了する。

      $ cd $CURP_HOME/docs/_api/
      $ open api.pdf

   ただし、pdfを生成する場合にはさらにLaTeXが必要なので、
   インストールされていることを確認しておくこと。

7. 最新のcurpを取得する。

   curpのディレクトリに移動して、次のようにコマンドを打つ::

      $ hg pull # BitBucketサーバからソースの最新版の履歴を取得。
      (パスワードを聞かれる)
      $ hg update # 履歴の最新版のファイルにアップデートする。

   もし, SSHの公開鍵を自分のアカウントに追加している場合にはパスワードを
   聞かれることはない．

以上の作業によって、最新のCurpを使うことが出来るようになる。
4.と5.で行った、Python環境の構築とfortranライブラリのビルドは基本的には
必要無い。
ただし、fortranライブラリが変更されている場合には、再ビルドを行うこと。

Try to run
===========

現時点では次のように実行する::

   $ $CURP_HOME/buildout/bin/curp-python $CURP_HOME/src/curp.py <input_file> > log

いくつかの例が$CURP_HOME/exampleにあるので、run.shを実行して試してみると良い。
