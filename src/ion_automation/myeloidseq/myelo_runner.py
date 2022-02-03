from myelo_worker import myeloseq

myelo_runner = myeloseq("/home/ionadmin/ion_config.conf")
myelo_runner.workbook = "/mnt/Z_drive/Molecular/IonTorrent/myeloseqer_test/worksheets/22-MGMQ1.xlsm"
myelo_runner.start()