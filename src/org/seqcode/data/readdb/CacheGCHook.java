package org.seqcode.data.readdb;

import java.util.logging.*;

/**
 * wait around and run garbage collection
 */
public class CacheGCHook implements Runnable {
	private Logger logger;
    private int gcFrequency;
    public CacheGCHook(Logger l, int gcFreq) {
        logger = l;
        gcFrequency = gcFreq;
    }

    public void run() {
        while (true) {
            if (LRUCache.removed() > gcFrequency) {
                logger.log(Level.INFO,"running GC");
                LRUCache.resetRemoved();
                System.gc();
                System.runFinalization();
            } 
            try {
                Thread.sleep(500);
            } catch (InterruptedException e) {

            }
        }
    }

}