import smtplib
import os
import csv
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
import cPickle


def send_email(my_obj, f_name):
    # create message
    msg = MIMEMultipart('alternative')

    # generate a pickled object
    pik_ob = cPickle.dumps(my_obj)

    # attach the object file
    filename = f_name
    # f = file(pik_ob)
    attachment = MIMEText(pik_ob)
    attachment.add_header('Content-Disposition', 'attachment', filename=filename)
    msg.attach(attachment)

    # add content to log message
    msg['Subject'] = "This is the subject"
    msg['From'] = 'pythonsend123@gmail.com'
    body = """This is the body of the message"""
    content = MIMEText(body, 'plain')
    msg.attach(content)

    # email the message
    server = smtplib.SMTP('smtp.gmail.com:587')
    server.starttls()
    server.login('pythonsend123@gmail.com', 'pythonsend')
    server.sendmail('pythonsend123@gmail.com', 'ewan.pinnington@gmail.com', msg.as_string())
    server.quit()