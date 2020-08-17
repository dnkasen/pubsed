#!/usr/bin/python
# script to run nightly regression test and
# put results in dropbox and email people
# 
# To start dropbox daemon try
# ~/.dropbox-dist/dropboxd
# To check status
# top | grep dropbox



import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.mime.text import MIMEText
from email.utils import COMMASPACE, formatdate
from email import encoders
import os.path as op
import os
import time
import subprocess


machine = "quest"
nprocs  = "8"
#send_to = ['daniel.kasen@gmail.com']
send_to = ['daniel.kasen@gmail.com','nroth@astro.umd.edu','hklion@berkeley.edu','abigail@berkeley.edu','david_khatami@berkeley.edu','kenshen@astro.berkeley.edu','dewillcox@lbl.gov','sricher@ncsu.edu','bennytth@gmail.com','siva.darbha@berkeley.edu']
result_dir = "/home/kasen/Dropbox/Sedona_Test_Results/"

def main():

	sedona_dir = os.environ['SEDONA_HOME']
	src_dir  = sedona_dir + "/src"
	test_dir = sedona_dir + "/tests"
	files = []
	compile_failed = 0

	testfile_pdf = 'test_results_' + time.strftime("%m-%d-%y") + '.pdf'
	testfile_txt = 'test_results_' + time.strftime("%m-%d-%y") + '.txt'
	message = "Sedona Nightly Test\n----------------\n\n"
        message += "To see test reults, go here:\n"
        message += "https://www.dropbox.com/sh/oaxypp1kov9babb/AAAwuduf0QRzhOytISQiistga?dl=0 \n\n"



	# pull the latest version
	os.chdir(sedona_dir)
	cmd = "git pull origin development"
	message += "running: "  + cmd + "\n"
	message += "---------------------------------------\n"
	try:
		result = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
		message = message + result + "\n"
	except subprocess.CalledProcessError as e:
		message += "\GitHub pull failed. Output follows:\n"
		message += "-------------------------------------\n"
		message += e.output + "\n"
		compile_failed = 1
		subject = "Sedona Test: GitHub pull FAILED"


	# make the latest version
	os.chdir(src_dir)
	cmd  = "./install.sh clean"
	result = subprocess.check_output(cmd, shell=True)
	cmd  = "./install.sh " + machine
	try:
		subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
	except subprocess.CalledProcessError as e:
		message += "\nCompilation failed. Output follows:\n"
		message += "-------------------------------------\n"
		message += e.output + "\n"
		compile_failed = 2
		subject = "Sedona Test: Compilation FAILED"


	# if compilation didn't fail, run tests
	if (not compile_failed):
		os.chdir(test_dir)
		cmd = "python run_test_suite.py -n " + nprocs

		try:
			result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
			message = message + "\n\n" + result
#			files = [testfile_pdf,testfile_txt]
#			files = [testfile_pdf]
                        files = []

			# determine if passed/failed
			if ('FAILED' in message):
				subject = "Sedona Test: FAILED"
			else:
				subject = "Sedona Test: PASSED"
		except subprocess.CalledProcessError as e:
			subject = "Sedona Test: Runs Crashed"
			message += "\nRun script failed: Output follows:\n"
			message += "-------------------------------------\n"
			message += e.output + "\n"

	# send out the results
	email_results(send_to,message,subject,files)

	# copy over results
	os.system("mv test_results_* " + result_dir)


def email_results(send_to,message,subject,files):

	user = 'sedona.transport@gmail.com'

	msg = MIMEMultipart()
	msg['Subject'] = subject
	msg['From'] = user
	msg['To'] = COMMASPACE.join(send_to)
	msg['Date'] = formatdate(localtime=True)
	msg.attach(MIMEText(message))


	for path in files:
		part = MIMEBase('application', "octet-stream")
		with open(path, 'rb') as file:
			part.set_payload(file.read())
		encoders.encode_base64(part)
		part.add_header('Content-Disposition','attachment; filename="{}"'.format(op.basename(path)))
		msg.attach(part)


	server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
	server.ehlo()
	server.login(user, PUT_PASSWORD_HERE)
	server.sendmail(user, send_to, msg.as_string())
	server.close()

main()
