tar -cf - GAP_0AA_600AA_AF3_SCREEN | pigz -9 -p $(nproc) > GAP_0AA_600AA_AF3_SCREEN.tar.gz
tar -cf - GAP_600AA_1000AA_AF3_SCREEN | pigz -9 -p $(nproc) > GAP_600AA_1000AA_AF3_SCREEN.tar.gz
tar -cf - GAP_1000AA_1500AA_AF3_SCREEN | pigz -9 -p $(nproc) > GAP_1000AA_1500AA_AF3_SCREEN.tar.gz
tar -cf - GAP_1500AA_2000AA_AF3_SCREEN | pigz -9 -p $(nproc) > GAP_1500AA_2000AA_AF3_SCREEN.tar.gz

tar -cf - GEF_0AA_600AA_AF3_SCREEN | pigz -9 -p $(nproc) > GEF_0AA_600AA_AF3_SCREEN.tar.gz
tar -cf - GEF_600AA_1000AA_AF3_SCREEN | pigz -9 -p $(nproc) > GEF_600AA_1000AA_AF3_SCREEN.tar.gz
tar -cf - GEF_1000AA_1500AA_AF3_SCREEN | pigz -9 -p $(nproc) > GEF_1000AA_1500AA_AF3_SCREEN.tar.gz
tar -cf - GEF_1500AA_2000AA_AF3_SCREEN | pigz -9 -p $(nproc) > GEF_1500AA_2000AA_AF3_SCREEN.tar.gz
