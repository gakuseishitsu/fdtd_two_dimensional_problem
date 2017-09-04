#-------------------------------------------------------------------------------
# gnuplot�̐ݒ�
#-------------------------------------------------------------------------------
reset
set nokey                 # �}��̔�\��
set xrange [0:0.6]       # x�������͈̔͂̐ݒ�
set yrange [0:0.6]         # y�������͈̔͂̐ݒ�
set zrange [-1500:1500]         # z�������͈̔͂̐ݒ�

set term gif animate delay 15      # �o�͂�gif�A�j���ɐݒ�
set output "fdtd_2.gif"  # �o�̓t�@�C�����̐ݒ�

#set contour
#set view 0,0
#unset surface
#set cntrparam levels 10

#set object 1 rect from 1.0, -0.4 to 1.1, 1 back linewidth 0 fillcolor rgb "#dcdcdc" fill solid 0.5
#set object 2 rect from 0, -0.4 to 0.008, 1 back linewidth 0 fillcolor rgb "#dcdcdc" fill solid 0.5
#set object 3 rect from 2, -0.4 to 2.016, 1 back linewidth 0 fillcolor rgb "#dcdcdc" fill solid 0.5

#-------------------------------------------------------------------------------
# �ϐ��̐ݒ�
#-------------------------------------------------------------------------------
n0 = 1        # ���[�v�ϐ��̏����l
n1 = 70     # ���[�v�ϐ��̍ő�l
dn = 1  # ���[�v�ϐ��̑����Ԋu

#-------------------------------------------------------------------------------
# ���[�v�̊J�n
#-------------------------------------------------------------------------------
load "gnuplot_fdtd_2.plt" 