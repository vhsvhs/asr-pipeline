"""
    Tools for interfacing with Amazon Web Services S3 storage.
    These methods are optional, and can be enabled using --enable_aws
    at the command line when running the ASR pipeline.
"""

import boto
import boto.ec2
import boto.s3
import boto.sqs
from boto.sqs.message import Message
from boto.s3.connection import S3Connection
from boto.s3.connection import Location
import sys, time

S3LOCATION = Location.USWest
ZONE = "us-west-1"

def aws_update_status(message, S3_BUCKET, S3_KEYBASE):
    if S3_BUCKET == None:
        print "\n. Error, you haven't defined an S3 bucket."
        exit()
    if S3_KEYBASE == None:
        print "\n. Error, you haven't defined an S3 key base."
        exit()
                
    """Update the status field in AWS S3 for this job"""
    #s3 = boto.connect_s3()
    s3 = S3Connection()
    
    bucket = s3.lookup(S3_BUCKET)
    if bucket == None:
        bucket = s3.create_bucket(S3_BUCKET, location=S3LOCATION)
        bucket.set_acl('public-read')
    
    STATUS_KEY = S3_KEYBASE + "/status"
    key = bucket.get_key(STATUS_KEY)
    if key == None:
        key = bucket.new_key(STATUS_KEY)

        #key = bucket.get_key(STATUS_KEY) 
    if key == None:
        print "\n. Error 39 - the key is None"
        exit()   
    key.set_contents_from_string(message)
    key.set_acl('public-read')
    
    print "\n. S3 Status Update:", key.get_contents_as_string()

def aws_checkpoint(checkpoint, S3_BUCKET, S3_KEYBASE):
    if S3_BUCKET == None:
        print "\n. Error, you haven't defined an S3 bucket."
        exit()
    if S3_KEYBASE == None:
        print "\n. Error, you haven't defined an S3 key base."
        exit()
                
    """Update the status field in AWS S3 for this job"""
    #s3 = boto.connect_s3()
    s3 = S3Connection()
    
    bucket = s3.lookup(S3_BUCKET)
    if bucket == None:
        bucket = s3.create_bucket(S3_BUCKET, location=S3LOCATION)
        bucket.set_acl('public-read')
    
    CHECKPOINT_KEY = S3_KEYBASE + "/checkpoint"
    key = bucket.get_key(CHECKPOINT_KEY)
    if key == None:
        key = bucket.new_key(CHECKPOINT_KEY)
        #key = bucket.get_key(STATUS_KEY) 
    if key == None:
        print "\n. Error 67 - the key is None"
        exit()   
    key.set_contents_from_string(checkpoint.__str__())
    key.set_acl('public-read')    
    print "\n. S3 Checkpoint:", key.get_contents_as_string()

def push_database_to_s3(dbpath, S3_BUCKET, S3_KEYBASE):
    """Pushes the startup files for a job to S3.
        jobid is the ID of the Job object,
        filepath is the path on this server to the file."""
    #s3 = S3Connection()
    print "\n. Pushing the database to S3"
    
    s3 = S3Connection()
    
    #s3 = boto.connect_s3()
    #print "46: connected to", s3.Location
    
#     allBuckets = s3.get_all_buckets()
#     for bucket in allBuckets:
#         print(str(bucket.name))
#         allKeys = bucket.get_all_keys()
#         for key in allKeys:
#             print "\t", key.name
    
    bucket = s3.lookup(S3_BUCKET)
    if bucket == None:
        bucket = s3.create_bucket(S3_BUCKET, location=S3LOCATION)
        bucket.set_acl('public-read')
    
    
    SQLDB_KEY = S3_KEYBASE + "/sqldb"
    key = bucket.get_key(SQLDB_KEY)
    if key == None:
        key = bucket.new_key(SQLDB_KEY)
    key.set_contents_from_filename(dbpath)
    key.set_acl('public-read')    
    
    
def sqs_stop(jobid, attempts=0):
    conn = boto.sqs.connect_to_region(ZONE)
    queue = conn.get_queue("phylobot-jobs")
    if queue == None:
        queue = conn.create_queue("phylobot-jobs")
    m = Message()
    m.set_body('stop ' + jobid.__str__() + " " + attempts.__str__())
    queue.write(m)  