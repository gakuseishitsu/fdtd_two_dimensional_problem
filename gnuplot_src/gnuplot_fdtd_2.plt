# gnuplot_fdtd_2.plt

#-------------------------------------------------------------------------------
# ���[�v����
#-------------------------------------------------------------------------------
if(exist("n")==0 || n<0) n = n0  # ���[�v�ϐ��̏�����

#-------------------------------------------------------------------------------
# �v���b�g
#-------------------------------------------------------------------------------
plot "fdtd_data.dat"  index n using 1:2 with lines  lc rgb "#000000" lw 2 # n�Ԗڂ̃f�[�^�̃v���b�g

#-------------------------------------------------------------------------------
# ���[�v����
#-------------------------------------------------------------------------------
n = n + dn            # ���[�v�ϐ��̑���
if ( n < n1 ) reread  # ���[�v�̕]��
undefine n            # ���[�v�ϐ��̍폜