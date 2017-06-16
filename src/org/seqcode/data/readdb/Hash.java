package org.seqcode.data.readdb;
import java.security.*;

public class Hash {
	public static String getHash(String txt, String hashType) {
        try {
            MessageDigest md = MessageDigest.getInstance(hashType);
            byte[] array = md.digest(txt.getBytes());
            StringBuffer sb = new StringBuffer();
            for (int i = 0; i < array.length; ++i) {
            	sb.append(Integer.toHexString((array[i] & 0xFF) | 0x100).substring(1,3));
            }
            return sb.toString();
        } catch (NoSuchAlgorithmException e) {
            //error action
        }
        return null;
    }

    public static String md5(String txt) {
        return getHash(txt, "MD5");
    }
}