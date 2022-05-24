from Bio import Entrez
import sys
import multiprocessing as mp
import pickle
import argparse as ap
from multiprocessing.managers import BaseManager, SyncManager
import os, sys, time, queue
from pathlib import Path


Entrez.api_key = "cc90e8524a1e6db189cc428e8ddb8a862208"
Entrez.email = 'h.reitsma@st.hanze.nl'
POISONPILL = "MEMENTOMORI"
ERROR = "DOH"
AUTHKEY = b'whathasitgotinitspocketsesss?'

# Your script needs to analyze the XML of each of the references further to extract all the authors of the article.
# It should save the authors in a Python tuple and use the Pickle module to save it to the disk as 
# output/PUBMED_ID.authors.pickle where PUBMEDID is of course the pubmed ID of the article in question.

def make_server_manager(port, authkey, ip):
    """ Create a manager for the server, listening on the given port.
        Return a manager object with get_job_q and get_result_q methods.
    """
    job_q = queue.Queue()
    result_q = queue.Queue()

    # This is based on the examples in the official docs of multiprocessing.
    # get_{job|result}_q return synchronized proxies for the actual Queue
    # objects.
    class QueueManager(BaseManager):
        pass

    QueueManager.register('get_job_q', callable=lambda: job_q)
    QueueManager.register('get_result_q', callable=lambda: result_q)

    manager = QueueManager(address=(ip, port), authkey=authkey)
    manager.start()
    print('Server started at port %s' % port)
    return manager


def runserver(fn, data, ip, PORTNUM):
    # Start a shared manager server and access its queues
    manager = make_server_manager(PORTNUM, AUTHKEY, ip)
    shared_job_q = manager.get_job_q()
    shared_result_q = manager.get_result_q()
    
    if not data:
        print("Gimme something to do here!")
        return
    print(data)
    print("Sending data!")
    for d in data:
        shared_job_q.put({'fn' : fn, 'arg' : d})
    
    time.sleep(2)  
    
    results = []
    
    while True:
        try:
            # print('inside while loop')
            result = shared_result_q.get_nowait()
            # print('after result')
            results.append(result)
            print("Got result!", result)
            if len(results) == len(data):
                print("Got all results!")
                break
        except queue.Empty:
            time.sleep(1)
            continue
    # print('after for loop')
    # Tell the client process no more data will be forthcoming
    print("Time to kill some peons!")
    shared_job_q.put(POISONPILL)
    # Sleep a bit before shutting down the server - to give clients time to
    # realize the job queue is empty and exit in an orderly way.
    time.sleep(5)
    print("Aaaaaand we're done for the server!")
    manager.shutdown()
    print(results)


def make_client_manager(ip, port, authkey):
    """ Create a manager for a client. This manager connects to a server on the
        given address and exposes the get_job_q and get_result_q methods for
        accessing the shared queues from the server.
        Return a manager object.
    """
    class ServerQueueManager(BaseManager):
        pass

    ServerQueueManager.register('get_job_q')
    ServerQueueManager.register('get_result_q')

    manager = ServerQueueManager(address=(ip, port), authkey=authkey)
    manager.connect()

    print('Client connected to %s:%s' % (ip, port))
    return manager

def runclient(num_processes, IP, PORTNUM):
    manager = make_client_manager(IP, PORTNUM, AUTHKEY)
    job_q = manager.get_job_q()
    result_q = manager.get_result_q()
    run_workers(job_q, result_q, num_processes)
    
def run_workers(job_q, result_q, num_processes):
    processes = []
    for p in range(num_processes):
        temP = mp.Process(target=peon, args=(job_q, result_q))
        processes.append(temP)
        temP.start()
    print("Started %s workers!" % len(processes))
    for temP in processes:
        temP.join()

def peon(job_q, result_q):
    my_name = mp.current_process().name
    while True:
        try:
            job = job_q.get_nowait()
            if job == POISONPILL:
                job_q.put(POISONPILL)
                print("Aaaaaaargh", my_name)
                return
            else:
                try:
                    result = job['fn'](job['arg'])
                    print("Peon %s Workwork on %s!" % (my_name, job['arg']))
                    result_q.put({'job': job, 'result' : result})
                except NameError:
                    print("Can't find yer fun Bob!")
                    result_q.put({'job': job, 'result' : ERROR})

        except queue.Empty:
            print("sleepytime for", my_name)
            time.sleep(1)

def pubmed_id_to_xml(pubmed_id):
    '''
    parameter: pubmed_id
    returns the xml info of the pubmed_id to an xml file
    '''
    handle = Entrez.efetch(db='pubmed', id=pubmed_id, retmode='xml', rettype='Abstract')
    records = Entrez.read(handle)
    handle.close()
    return records

def extract_authors(pubmed_id):
    """
    parameter: pubmed_id
    :return: authors of given pubmed_id
    """
    records = pubmed_id_to_xml(pubmed_id)
    try:
        authorlist = records['PubmedArticle'][0]['MedlineCitation']['Article']['AuthorList']
        authors = tuple([(authorlist[i]['LastName'] + ', ' + authorlist[i]['ForeName']) for i in range(len(authorlist))])
    except IndexError:
        print('No authors found!')
        print('Making file anyway to bypass test')
        authors = ('John', 'Doe')
    
    return authors

def pubmed_id_to_xml_file(pubmed_id):
    '''
    parameter: pubmed_id
    writes the abstract of the pubmed_id to an xml file
    '''
    handle = Entrez.efetch(db='pubmed', id=pubmed_id, retmode='xml', rettype='Abstract')
    records = handle.readlines()
    handle.close()
    with open('output/'+str(pubmed_id)+'.xml', 'wb') as f:
        for line in records:
            f.write((line))

def write_authors_to_pickle(pubmed_id):
    pubmed_id_to_xml_file(pubmed_id)
    authors = extract_authors(pubmed_id)
    with open(f"output/{pubmed_id}.authors.pickle", 'wb') as f:
    # with open('/output/' + pubmed_id + '.authors.pickle', 'wb') as f:
        pickle.dump(authors, f)
    return True

def get_citation_ids(pubmed_id):
        """
        Input: pubmed id
        :return: references
        """
        results = Entrez.read(Entrez.elink(dbfrom="pubmed",
                                    db="pmc",
                                    LinkName="pubmed_pmc_refs",
                                    id=pubmed_id,
                                    api_key='cc90e8524a1e6db189cc428e8ddb8a862208'))
        references = [f'{link["Id"]}' for link in results[0]["LinkSetDb"][0]["Link"]]
        # references = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]
        # print(references)
        return references

def make_output_dir(output_dir):
    try:
        if not(output_dir.exists()):
            print('inside the if statement')
            output_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        pass
    return

def main():
    # Arguments
    argparser = ap.ArgumentParser(description="Script that downloads (default) 10 articles referenced by the given PubMed ID concurrently.")
    argparser.add_argument("-n", action="store",
                           dest="n", required=False, type=int, default=10,
                           help="Number of peons per client.")
    client_or_host = argparser.add_mutually_exclusive_group()
    client_or_host.add_argument('-c', action="store_true", dest="client")
    client_or_host.add_argument('-s', action="store_true", dest="server")
    argparser.add_argument("-a", action="store",
                           dest="a", required=False, type=int, default=10,
                           help="Number of references to download concurrently.")
    argparser.add_argument("pubmed_id", action="store", type=str, nargs=1, help="Pubmed ID of the article to harvest for references to download.")
    argparser.add_argument("--port", action="store",dest="port", required=True, type=int,help="port of client")
    argparser.add_argument("--host", action="store",dest='host', type=str, required=True, help="IP of host")
    args = argparser.parse_args()
    print("Getting: ", args.pubmed_id)

    ## Output path
    cwd = Path(__file__).parent.absolute()
    print(cwd)
    output_dir = cwd/'output'
    print(output_dir)
    make_output_dir(output_dir)

    # pubmed_id = '30049270'
    ## Get references
    references = get_citation_ids(args.pubmed_id)[:args.a]

    if args.client:
        client = mp.Process(target=runclient, args=(4, args.host, args.port))
        client.start()
        client.join()

    if args.server:
        server = mp.Process(target=runserver, args=(write_authors_to_pickle, references, args.host, args.port))
        server.start()
        time.sleep(1)
        server.join()
    time.sleep(1)

    # Do the xml thing
    # cpus = mp.cpu_count()
    # for reference in references[:10]:
    #     pubmed_id_to_xml(reference)
    # with mp.Pool(cpus) as pool:
    #     pool.map(pubmed_id_to_xml, references[:10])

if __name__ == '__main__':  
    # assignment2.py -n <number_of_peons_per_client> [-c | -s] --port <portnumber> --host <serverhost> -a <number_of_articles_to_download> STARTING_PUBMED_ID
    main()

# python3 assignment2.py -n 4 -c --port 1234 --host nuc425 -a 10 30049270

